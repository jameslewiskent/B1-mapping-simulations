function [simulation_results,settings] = epg_dream(settings)
% Function to simulate image trains for DREAM B1 mapping scheme
% Controls looping structure over different flip angles, T1s etc.
%
%   INFO:
%       Dual refocussing echo acquisition mode (DREAM) is an ultra-fast B1
%       mapping scheme using the images formed from a STEAM preparation
%       sequence followed by a fast image train. The ratio of the FID and
%       STE signals gives the actual flip angle from:
%
%       Actual FA = atan(sqrt((2*IT1/IT2)))
%           Eq. [5] from Nehrke, K. et al. MRM 68:1517-1526 (2012).
%           https://doi.org/10.1002/mrm.24158.
%
% J. Kent. 2023. Using B.Hargreaves EPG and Bloch Functions.

% ----------------------- Synthetic Data Sim ----------------------- %
%IT1 = zeros(settings.Train_Size,size(settings.slice_T1,1),size(settings.slice_T1,2),size(settings.Tx_FA_map,3)); % T1 Dimension set to 2 to prevent collapse
%IT2 = zeros(settings.Train_Size,size(settings.slice_T1,1),size(settings.slice_T1,2),size(settings.Tx_FA_map,3)); % T1 Dimension set to 2 to prevent collapse
if settings.UseSyntheticData == 1
    % Go through image indexes
    %     for indi = 1:size(slice_T1,1)
    %         for indj = 1:size(slice_T1,2)
    %             if slice_T1(indi,indj) ~= 0 % Only run Bloch simulation for heart & blood (non-zero T1s), ignore rest of anatomy
    %                 % Gets T1,T2 and PP_FA for voxel
    %                 T1 = slice_T1(indi,indj);
    %                 T2 = slice_T2(indi,indj);
    %                 PP_FA = Tx_FA_map(indi,indj).*pi/180; % Pre-pulse Flip Angles in Radians
    %                % for ratio_n = 1:length(ratios)
    %                     IT_FA = (PP_FA./ratios(1)).*pi/180;  % Image Train Flip Angles in Radians
    %                     % for now only uses first ratio in array
    %
    %                     for TR_n = 1:length(TRs)
    %                         TR = TRs(TR_n);
    %
    %                         [IT1(indi,indj,TR_n,1,:),IT2(indi,indj,TR_n,1,:)] = DREAM_Seq(PP_FA,IT_FA,T1,T2,Td,Ts,phi1,phi2,kg,Diff_co,Velocity,Angle,prep_spoils,noadd);
    %
    %                     end
    %                 %end
    %             end
    %         end
    %     end
    %     PP_FAs = 'none';
    %     IT_FAs = 'none';
else
    IT1 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3),length(settings.Dynamic_Range),length(settings.B0_Range_Hz),length(settings.T1s));
    IT2 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3),length(settings.Dynamic_Range),length(settings.B0_Range_Hz),length(settings.T1s));
    Total_Energy = zeros(1,length(settings.Dynamic_Range));
    for Dynamic_Range_n = 1:length(settings.Dynamic_Range)
        for T1_n = 1:length(settings.T1s)
            T1 = settings.T1s(T1_n);
            for B0_n = 1:length(settings.B0_Range_Hz)
                for Flow_n = 1:length(settings.Velocities)
                    Velocity = settings.Velocities(Flow_n);
                    Angle = settings.Angles(Flow_n);
                    for Diff_n = 1:length(settings.Diff_coeffs)
                        Diff_co = settings.Diff_coeffs(Diff_n);
                        % Simulate B0 effect on FA
                        IT_FAs = Simulate_RF_Pulse(settings.Hz_per_Volt*settings.Dynamic_Range(Dynamic_Range_n)*settings.IT_RF_Pulse,settings.IT_RF_Time,settings.B0_Range_Hz(B0_n)+settings.Slice_Shifts,T1,settings.T2,settings.Gamma); % Output FA in radians
                        PP_FAs = Simulate_RF_Pulse(settings.Hz_per_Volt*settings.Dynamic_Range(Dynamic_Range_n)*settings.PP_RF_Pulse,settings.PP_RF_Time,settings.B0_Range_Hz(B0_n)+settings.PP_Shifts,T1,settings.T2,settings.Gamma); % Output FA in radians
                        settings = Check_Mag_Track_Flag(Dynamic_Range_n,T1_n,B0_n,Flow_n,Diff_n,settings);
                        [IT1(:,:,:,Dynamic_Range_n,B0_n,T1_n,Flow_n,Diff_n),IT2(:,:,:,Dynamic_Range_n,B0_n,T1_n,Flow_n,Diff_n),Cumulative_Time,N_imaging_RF,N_prep_RF,Mag_Track] = DREAM_Seq(T1,IT_FAs,PP_FAs,Velocity,Angle,Diff_co,settings);
                        if any(settings.Mag_Track_Flags == 1)
                            simulation_results.Mag_Track{settings.Mag_Track_Flags == 1} = Mag_Track;
                        end
                    end
                end
            end
        end
        Imaging_RF_Energy = Calc_RF_Energy(settings.Dynamic_Range(Dynamic_Range_n)*settings.IT_RF_Pulse,settings.IT_RF_Time); % Calculate Imaging RF pulse energy
        Preparation_RF_Energy = Calc_RF_Energy(settings.Dynamic_Range(Dynamic_Range_n)*settings.PP_RF_Pulse,settings.PP_RF_Time); % Calculate Preparation RF pulse energy
        Total_Energy(Dynamic_Range_n) = N_imaging_RF*Imaging_RF_Energy + N_prep_RF*Preparation_RF_Energy;
    end
end
Average_10s_Power = 10*Total_Energy./Cumulative_Time;

% Sequence Function
    function [Train1,Train2,Cumulative_Time,N_imaging_RF,N_prep_RF,Mag_Track] = DREAM_Seq(T1,IT_FAs,PP_FAs,Velocity,Angle,Diff_co,settings)
        % Function runs inside relevent epg_SCHEME function and simulates a single set of image trains for specified parameters.
        
        P = [0 0 1]';% Z0 = 1 (Equilibrium)
        Mag_Track(:,1) = [P(3,1);0];
        Train1 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3));
        Train2 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3));
        
        Cumulative_Time = 0;
        N_prep_RF = 0;
        N_imaging_RF = 0;
        
        
        if strcmp(settings.MSor3D,'3D')
            %%% START OF 3D SEQUENCE %%%
            for Mode_n = 1:size(settings.Tx_FA_map,3) % repeat for different modes of mTx array
                
                % STEAM PREPARATION
                N_prep_RF = N_prep_RF + 1;
                [P,~,Mag_Track] = epg_rf(P,PP_FAs(settings.Centre_Slice),settings.RF_Phase(Mode_n),settings.PP_RF_Time,Mag_Track,settings); % 1st STEAM RF alpha pulse
                
                time1 = settings.Ts/2;
                Cumulative_Time = Cumulative_Time + time1;
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,time1,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation
                
                if strcmp(settings.Echo_Order,'STEFirst')
                    P = epg_mgrad(P); % Adds negative "unit" Gm2 gradient
                elseif strcmp(settings.Echo_Order,'FIDFirst')
                    P = epg_grad(P); % Adds positive "unit" Gm2 gradient
                end
                
                Cumulative_Time = Cumulative_Time + time1;
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,time1,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation
                
                N_prep_RF = N_prep_RF + 1;
                [P,~,Mag_Track] = epg_rf(P,PP_FAs(settings.Centre_Slice),settings.RF_Phase(Mode_n),settings.PP_RF_Time,Mag_Track,settings); % 2nd STEAM RF alpha pulse
                
                time2 = settings.Td - settings.Ts;
                Cumulative_Time = Cumulative_Time + time2;
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,time2,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation
                
                for n=1:settings.prep_spoils
                    P = epg_grad(P,settings.noadd); % spoiling
                end
                % END OF STEAM PREPARATION
                
                % IMAGING PULSE TRAIN
                for PE2_n = 1:settings.Scan_Size(2)
                    Slice_n = settings.Slice_Order(PE2_n);
                    for IT_n = 1:settings.Scan_Size(1)
                        N_imaging_RF = N_imaging_RF + 1;
                        [P,~,Mag_Track] = epg_rf(P,IT_FAs(Slice_n),settings.RF_Phase(Mode_n)+settings.rf_spoiling_phases(N_imaging_RF),settings.IT_RF_Time,Mag_Track,settings); % Imaging RF
                        
                        % STE FIRST FID SECOND SCHEME (Fig 1.a Nehkre, K. 2012)
                        Cumulative_Time = Cumulative_Time + settings.TE1;
                        [P,~,~,Mag_Track] =  epg_grelax(P,T1,settings.T2,settings.TE1,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation time = TE_STE
                        
                        P = epg_grad(P); % Adds positive "unit" Gm gradient
                        
                        Train1(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*settings.rf_spoiling_phases(N_imaging_RF)); % Store Phase-Demodulated 'signal' ; %P(2,2); % record STE 'signal'
                        
                        Cumulative_Time = Cumulative_Time + (settings.TE2 - settings.TE1);
                        [P,~,~,Mag_Track] =  epg_grelax(P,T1,settings.T2,settings.TE2 - settings.TE1,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation time
                        
                        P = epg_mgrad(P); % Adds negative "unit" Gm gradient
                        
                        Train2(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*settings.rf_spoiling_phases(N_imaging_RF)); % Store Phase-Demodulated 'signal'; % record FID 'signal'
                        
                        for i = 1:10
                            P = epg_grad(P,settings.noadd); % Spoiling
                        end
                        
                        Cumulative_Time = Cumulative_Time + (settings.TR - settings.TE2);
                        [P,~,~,Mag_Track] =  epg_grelax(P,T1,settings.T2,settings.TR - settings.TE2,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation time = TR - TE_FID
                    end
                end
            end
            %%% END OF 3D SEQUENCE %%%
        elseif strcmp(settings.MSor3D,'2D')
            %%% START OF 2D SEQUENCE %%%
            for Mode_n = 1:size(settings.Tx_FA_map,3) % repeat for different modes of mTx array
                for PE2_n = 1:settings.Scan_Size(2)
                    Slice_n = settings.Slice_Order(PE2_n);
                    % STEAM PREPARATION
                    N_prep_RF = N_prep_RF + 1;
                    [P,~,Mag_Track] = epg_rf(P,PP_FAs(Slice_n),settings.RF_Phase(Mode_n),settings.PP_RF_Time,Mag_Track,settings); % 1st STEAM RF alpha pulse
                    
                    time1 = settings.Ts/2;
                    Cumulative_Time = Cumulative_Time + time1;
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,time1,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation
                    
                    if strcmp(settings.Echo_Order,'STEFirst')
                        P = epg_mgrad(P); % Adds negative "unit" Gm2 gradient
                    elseif strcmp(settings.Echo_Order,'FIDFirst')
                        P = epg_grad(P); % Adds positive "unit" Gm2 gradient
                    end
                    
                    Cumulative_Time = Cumulative_Time + time1;
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,time1,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation
                    
                    N_prep_RF = N_prep_RF + 1;
                    [P,~,Mag_Track] = epg_rf(P,PP_FAs(Slice_n),settings.RF_Phase(Mode_n),settings.PP_RF_Time,Mag_Track,settings); % 2nd STEAM RF alpha pulse
                    
                    time2 = settings.Td - settings.Ts;
                    Cumulative_Time = Cumulative_Time + time2;
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,time2,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation
                    
                    for n=1:settings.prep_spoils
                        P = epg_grad(P,settings.noadd); % spoiling
                    end
                    % END OF STEAM PREPARATION
                    
                    % IMAGING PULSE TRAIN
                    for IT_n = 1:settings.Scan_Size(1)
                        N_imaging_RF = N_imaging_RF + 1;
                        [P,~,Mag_Track] = epg_rf(P,IT_FAs(Slice_n),settings.RF_Phase(Mode_n)+settings.rf_spoiling_phases(N_imaging_RF),settings.IT_RF_Time,Mag_Track,settings); % Imaging RF
                        
                        % STE FIRST FID SECOND SCHEME (Fig 1.a Nehkre, K. 2012)
                        Cumulative_Time = Cumulative_Time + settings.TE1;
                        [P,~,~,Mag_Track] =  epg_grelax(P,T1,settings.T2,settings.TE1,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation time = TE_STE
                        
                        P = epg_grad(P); % Adds positive "unit" Gm gradient
                        
                        Train1(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*settings.rf_spoiling_phases(N_imaging_RF)); % Store Phase-Demodulated 'signal' ; %P(2,2); % record STE 'signal'
                        
                        Cumulative_Time = Cumulative_Time + (settings.TE2 - settings.TE1);
                        [P,~,~,Mag_Track] =  epg_grelax(P,T1,settings.T2,settings.TE2 - settings.TE1,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation time
                        
                        P = epg_mgrad(P); % Adds negative "unit" Gm gradient
                        
                        Train2(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*settings.rf_spoiling_phases(N_imaging_RF)); % Store Phase-Demodulated 'signal'; % record FID 'signal'
                        
                        for i = 1:10
                            P = epg_grad(P,settings.noadd); % Spoiling
                        end
                        
                        Cumulative_Time = Cumulative_Time + (settings.TR - settings.TE2);
                        [P,~,~,Mag_Track] =  epg_grelax(P,T1,settings.T2,settings.TR - settings.TE2,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation time = TR - TE_FID
                    end
                    % Period of relaxation
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.Compulsory_Delay_Time,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation + spoiling
                    Cumulative_Time = Cumulative_Time + settings.Compulsory_Delay_Time;
                end
            end
            %%% END OF 2D SEQUENCE %%%
        end
    end

simulation_results.IT1 = IT1;
simulation_results.IT2 = IT2;
simulation_results.N_imaging_RF = N_imaging_RF;
simulation_results.N_prep_RF = N_prep_RF;
simulation_results.Cumulative_Time = Cumulative_Time;
simulation_results.Total_Energy = Total_Energy;
simulation_results.Average_10s_Power = Average_10s_Power;
end
