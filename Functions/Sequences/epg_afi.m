function [simulation_results,settings] = epg_afi(settings)
% Function to simulate image trains for AFI B1 mapping scheme.
% Controls looping structure over different flip angles, T1s etc.
%
%   INFO:
%       Actual Flip Angle Imaging (AFI) is a B1 mapping scheme that relies
%       on perfect spoiling of the transverse magnetisation using two
%       reptition times of TR1 and TR2 << T1.
%
%       Actual FA = acos( (r.*n - 1) ./ (n - r) )
%           Eq. [6] in Yarnykh VL. AFI. Magn Reson Med 2007;57:192-200.
%           https://doi.org/10.1002/mrm.21120.
%
%       where n = TR2/TR1 and r is the ratio of the FID signals from TR1
%       and TR2.
%
% J. Kent. 2023. Using B.Hargreaves EPG and Bloch Functions.

% ----------------------- Synthetic Data Sim ----------------------- %
%IT1 = zeros(settings.Train_Size,size(settings.slice_T1,1),size(settings.slice_T1,2),2,size(settings.Tx_FA_map,3)); % T1 Dimension set to 2 to prevent collapse
%IT2 = zeros(settings.Train_Size,size(settings.slice_T1,1),size(settings.slice_T1,2),2,size(settings.Tx_FA_map,3)); % T1 Dimension set to 2 to prevent collapse
%Total_Energy = zeros(1,length(settings.Dynamic_Range));
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
    %                     FA = (PP_FA./ratios(1)).*pi/180;  % Image Train Flip Angles in Radians
    %                     % for now only uses first ratio in array
    %
    %                         [IT1(indi,indj,1,1,:),IT2(indi,indj,1,1,:)] = AFI_Seq(PP_FA,T1,T2,TR1,TR2,phi1,kg,Diff_co,Velocity,Angle,num_spoils,man_spoil,Dummy_Scans);
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
            for B0_n = 1:length(settings.B0_Range_Hz)
                T1 = settings.T1s(T1_n);
                FAs = Simulate_RF_Pulse(settings.Hz_per_Volt*settings.Dynamic_Range(Dynamic_Range_n)*settings.RF_Pulse,settings.RF_Time,settings.B0_Range_Hz(B0_n)+settings.Slice_Shifts,T1,settings.T2,settings.Gamma); % Output FA in radians
                settings = Check_Mag_Track_Flag(Dynamic_Range_n,T1_n,B0_n,settings);
                [IT1(:,:,:,Dynamic_Range_n,B0_n,T1_n),IT2(:,:,:,Dynamic_Range_n,B0_n,T1_n),Cumulative_Time,N_imaging_RF,Mag_Track] = AFI_Seq(T1,FAs,settings);
                if any(settings.Mag_Track_Flags == 1)
                    simulation_results.Mag_Track{settings.Mag_Track_Flags == 1} = Mag_Track;
                end
            end
        end
        Imaging_RF_Energy = Calc_RF_Energy(settings.Dynamic_Range(Dynamic_Range_n)*settings.RF_Pulse,settings.RF_Time); % Calculate Imaging RF pulse energy
        Total_Energy(Dynamic_Range_n) = N_imaging_RF*Imaging_RF_Energy;
    end
end
Average_10s_Power = 10*Total_Energy./Cumulative_Time;

% Sequence Function
    function [Train1,Train2,Cumulative_Time,N_imaging_RF,Mag_Track] = AFI_Seq(T1,FAs,settings)
        % Function runs inside relevent epg_SCHEME function and simulates a single set of image trains for specified parameters.
        
        P = [0 0 1]'; % Z0 = 1 (Equilibrium)
        Mag_Track(:,1) = [P(3,1);0];
        rTE1 = settings.TR1 - settings.TE;
        rTE2 = settings.TR2 - settings.TE;
        Train1 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3));
        Train2 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3));
        
        Cumulative_Time = 0;
        N_imaging_RF = 0;
        
        %%% START OF SEQUENCE %%%
        for Mode_n = 1:size(settings.Tx_FA_map,3) % repeat for different modes of mTx array
            
            % Dummy repetitions to produce steady state
            for Dummy_n = 1:settings.Dummy_Scans
                N_imaging_RF = N_imaging_RF +1;
                [P,~,Mag_Track] = epg_rf(P,FAs(settings.Slice_Order(Dummy_n)),settings.RF_Phase(Mode_n)+settings.rf_spoiling_phases(N_imaging_RF),settings.RF_Time,Mag_Track,settings);
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.TR1,settings.kg,settings.Diff_co,0,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation no spoiling
                Cumulative_Time = Cumulative_Time + settings.TR1;
                for i = 1:settings.num_spoils
                    P = epg_grad(P); % spoiling
                end
                if settings.man_spoil ==1
                    P = [0 0 P(3,1)]'; % manual spoiling
                end
                
                N_imaging_RF = N_imaging_RF +1;
                [P,~,Mag_Track] = epg_rf(P,FAs(settings.Slice_Order(Dummy_n)),settings.RF_Phase(Mode_n)+settings.rf_spoiling_phases(N_imaging_RF),settings.RF_Time,Mag_Track,settings);
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.TR2,settings.kg,settings.Diff_co,0,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation no spoiling
                Cumulative_Time = Cumulative_Time + settings.TR2;
                for i = 1:settings.num_spoils
                    P = epg_grad(P); % spoiling
                end
                if settings.man_spoil ==1
                    P = [0 0 P(3,1)]'; % manual spoiling
                end
            end
            
            for PE2_n = 1:settings.Scan_Size(2)
                Slice_n = settings.Slice_Order(PE2_n);
                for TR_n = 1:settings.Scan_Size(1)
                    N_imaging_RF = N_imaging_RF +1;
                    [P,~,Mag_Track] = epg_rf(P,FAs(Slice_n),settings.RF_Phase(Mode_n)+settings.rf_spoiling_phases(N_imaging_RF),settings.RF_Time,Mag_Track,settings);
                    
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.TE,settings.kg,settings.Diff_co,0,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation no spoiling
                    
                    Train1(TR_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*settings.rf_spoiling_phases(N_imaging_RF)); % Store Phase-Demodulated 'signal' F+0
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,rTE1,settings.kg,settings.Diff_co,0,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation no spoiling
                    Cumulative_Time = Cumulative_Time + settings.TR1;
                    for i = 1:settings.num_spoils
                        P = epg_grad(P); % spoiling
                    end
                    if settings.man_spoil ==1
                        P = [0 0 P(3,1)]'; % manual spoiling
                    end
                    N_imaging_RF = N_imaging_RF +1;
                    [P,~,Mag_Track] = epg_rf(P,FAs(Slice_n),settings.RF_Phase(Mode_n)+settings.rf_spoiling_phases(N_imaging_RF),settings.RF_Time,Mag_Track,settings);
                    
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.TE,settings.kg,settings.Diff_co,0,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation no spoiling
                    
                    Train2(TR_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*settings.rf_spoiling_phases(N_imaging_RF)); % Store Phase-Demodulated 'signal' F+0
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,rTE2,settings.kg,settings.Diff_co,0,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation no spoiling
                    Cumulative_Time = Cumulative_Time + settings.TR2;
                    for i = 1:settings.num_spoils
                        P = epg_grad(P); % spoiling
                    end
                    if settings.man_spoil ==1
                        P = [0 0 P(3,1)]'; % manual spoiling
                    end
                end
            end
        end
    end

simulation_results.IT1 = IT1;
simulation_results.IT2 = IT2;
simulation_results.N_imaging_RF = N_imaging_RF;
simulation_results.Cumulative_Time = Cumulative_Time;
simulation_results.Total_Energy = Total_Energy;
simulation_results.Average_10s_Power = Average_10s_Power;
end
