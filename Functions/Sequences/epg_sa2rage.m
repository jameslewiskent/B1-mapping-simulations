function [simulation_results,settings] = epg_sa2rage(settings)
% Function to simulate image trains for SA2RAGE B1 mapping scheme
% Controls looping structure over different flip angles, T1s etc.
%
%   INFO:
%
% J. Kent. 2023. Using B.Hargreaves EPG and Bloch Functions.


% End of checks and command output, following is control of loops through T1/TR/Ratio values etc.
% ----------------------- Synthetic Data Sim ----------------------- %
%IT1 = zeros(settings.Train_Size,size(settings.slice_T1,1),size(settings.slice_T1,2),2,size(settings.Tx_FA_map,3)); % T1 Dimension set to 2 to prevent collapse
%IT2 = zeros(settings.Train_Size,size(settings.slice_T1,1),size(settings.slice_T1,2),2,size(settings.Tx_FA_map,3)); % T1 Dimension set to 2 to prevent collapse
%Total_Energy = zeros(1,length(settings.Dynamic_Range));
if settings.UseSyntheticData == 1
    % Go through image indexes
    %     for indi = 1:size(settings.slice_T1,1)
    %         for indj = 1:size(settings.slice_T1,2)
    %             if settings.slice_T1(indi,indj) ~= 0 % Only run Bloch simulation for heart & blood (non-zero T1s), ignore rest of anatomy
    %                 % Gets T1,T2 and PP_FA for voxel
    %                 settings.T1 = settings.slice_T1(indi,indj);
    %                 settings.T2 = settings.slice_T2(indi,indj);
    %                 if settings.MTx == 0
    %                     PP_FA = abs(settings.Tx_FA_map(indi,indj)).*pi/180; % Pre-pulse Flip Angles in Radians
    %                     IT_FA = (PP_FA./FAratio);  % Image Train Flip Angles in Radians
    %                 else
    %                     % Else if mTx, pass mode values of PP_FA and IT_FA, sequence EPG
    %                     % function chooses value based on individual mode maps
    %                     PP_FA = abs(settings.Tx_FA_map(indi,indj,:)).*pi/180;
    %                     IT_FA = (PP_FA./FAratio);
    %                 end
    %                 RF_Phase = angle(conj(settings.Tx_FA_map(indi,indj,:))); % RF Phase from synthetic data (adds to RF spoiling in sim) (radians)
    %                 [IT1(:,indi,indj,:),IT2(:,indi,indj,:),Cumulative_Time)] = SA2RAGE_Seq(settings);
    %             end
    %         end
    %
    %         disp(['Simulating Image Trains. Percent Complete: ',num2str(100*indi/size(slice_T1,1),'%.2f'),'%'])
    %     end
    
else
    IT1 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3),length(settings.Dynamic_Range),length(settings.B0_Range_Hz),length(settings.T1s));
    IT2 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3),length(settings.Dynamic_Range),length(settings.B0_Range_Hz),length(settings.T1s));
    Total_Energy = zeros(1,length(settings.Dynamic_Range));
    for Dynamic_Range_n = 1:length(settings.Dynamic_Range)
        for T1_n = 1:length(settings.T1s)
            for B0_n = 1:length(settings.B0_Range_Hz)
                T1 = settings.T1s(T1_n);
                IT1_FAs = Simulate_RF_Pulse(settings.Hz_per_Volt*settings.Dynamic_Range(Dynamic_Range_n)*settings.IT1_RF_Pulse,settings.IT_RF_Time,settings.B0_Range_Hz(B0_n)+settings.Slice_Shifts,T1,settings.T2,settings.Gamma); % Output FA in radians
                IT2_FAs = Simulate_RF_Pulse(settings.Hz_per_Volt*settings.Dynamic_Range(Dynamic_Range_n)*settings.IT2_RF_Pulse,settings.IT_RF_Time,settings.B0_Range_Hz(B0_n)+settings.Slice_Shifts,T1,settings.T2,settings.Gamma); % Output FA in radians
                PP_FAs = Simulate_RF_Pulse(settings.Hz_per_Volt*settings.Dynamic_Range(Dynamic_Range_n)*settings.PP_RF_Pulse,settings.PP_RF_Time,settings.B0_Range_Hz(B0_n)+settings.PP_Shifts,T1,settings.T2,settings.Gamma); % Output FA in radians
                settings = Check_Mag_Track_Flag(Dynamic_Range_n,T1_n,B0_n,settings);
                [IT1(:,:,:,Dynamic_Range_n,B0_n,T1_n),IT2(:,:,:,Dynamic_Range_n,B0_n,T1_n),Cumulative_Time,N_imaging_RF1,N_imaging_RF2,N_prep_RF,Mag_Track] = SA2RAGE_Seq(T1,IT1_FAs,IT2_FAs,PP_FAs,settings);
                if any(settings.Mag_Track_Flags == 1)
                    simulation_results.Mag_Track{settings.Mag_Track_Flags == 1} = Mag_Track;
                end
            end
        end
        Imaging_RF1_Energy = Calc_RF_Energy(settings.Dynamic_Range(Dynamic_Range_n)*settings.IT1_RF_Pulse,settings.IT_RF_Time); % Calculate Imaging RF pulse energy
        Imaging_RF2_Energy = Calc_RF_Energy(settings.Dynamic_Range(Dynamic_Range_n)*settings.IT2_RF_Pulse,settings.IT_RF_Time); % Calculate Imaging RF pulse energy
        Preparation_RF_Energy = Calc_RF_Energy(settings.Dynamic_Range(Dynamic_Range_n)*settings.PP_RF_Pulse,settings.PP_RF_Time); % Calculate Preparation RF pulse energy
        Total_Energy(Dynamic_Range_n) = N_imaging_RF1*Imaging_RF1_Energy + N_imaging_RF2*Imaging_RF2_Energy + N_prep_RF*Preparation_RF_Energy;
    end
end
Average_10s_Power = 10*Total_Energy./Cumulative_Time;

% Sequence Function
    function [Train1,Train2,Cumulative_Time,N_imaging_RF1,N_imaging_RF2,N_prep_RF,Mag_Track] = SA2RAGE_Seq(T1,IT1_FAs,IT2_FAs,PP_FAs,settings)
        % Function runs inside relevent epg_SCHEME function and simulates a single set of image trains for specified parameters.
        
        P = [0 0 1]'; % Start magnetisation in equilibrium
        Mag_Track(:,1) = [P(3,1);0];
        IT_rTE = settings.IT_TR - settings.IT_TE;
        Train1 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3));
        Train2 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3));
        
        Cumulative_Time = 0;
        N_prep_RF = 0;
        N_imaging_RF1 = 0;
        N_imaging_RF2 = 0;
        
        %%% START OF SEQUENCE %%%
        if settings.Dummy_Scans ~= 0
            for Dummy_n = 1:settings.Dummy_Scans
                
                if size(settings.Tx_FA_map,3) > 1
                    Mode_n = size(settings.Tx_FA_map,3) - (settings.Dummy_Scans - Dummy_n); % Set mode for dummy scans to last in array
                    % PP_FA =
                    % IT_FA =
                else
                    Mode_n = 1; % Set mode for dummy scans
                end
                
                % Pre-pulse for dummy scans
                N_prep_RF = N_prep_RF +1;
                [P,~,Mag_Track] = epg_rf(P,PP_FAs(settings.Slice_Order(Dummy_n)),settings.RF_Phase(Mode_n),settings.PP_RF_Time,Mag_Track,settings);
                Cumulative_Time = Cumulative_Time + settings.PP_RF_Time;
                
                for i = 1:settings.prep_spoils(N_prep_RF)
                    P = epg_grad(P); % spoiling
                end
                if settings.man_spoil ==1
                    P = [0 0 P(3,1)]'; % manual spoiling
                end
                
                % Dummy image train 1
                for IT_n = 1:settings.Segment_Sizes(1)
                    N_imaging_RF1 = N_imaging_RF1 +1; % increment rf counter
                    [P,~,Mag_Track] = epg_rf(P,IT1_FAs(settings.Slice_Order(Dummy_n)),settings.RF_Phase(Mode_n)+settings.rf_spoiling_phases(N_imaging_RF1+N_imaging_RF2),settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                    
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TR,settings.kg,settings.Diff_co,0,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation no spoiling
                    
                    % Dummy Scans not stored
                    for i = 1:settings.train_spoils(IT_n)
                        P = epg_grad(P); % spoiling
                    end
                    if settings.man_spoil == 1
                        P = [0 0 P(3,1)]'; % manual spoiling
                    end
                    Cumulative_Time = Cumulative_Time + settings.IT_TR;
                end
                
                % Relaxation between dummy trains
                T_relax = settings.TD2 - settings.PP_RF_Time - 1.5*settings.Segment_Sizes(1).*settings.IT_TR; % Time for Relaxation between two trains(s)
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,T_relax,settings.kg,settings.Diff_co,1,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation + spoiling
                Cumulative_Time = Cumulative_Time + T_relax;
                
                % Dummy image train 2
                for IT_n = 1:settings.Segment_Sizes(1)
                    N_imaging_RF2 = N_imaging_RF2 +1; % increment rf counter
                    [P,~,Mag_Track] = epg_rf(P,IT2_FAs(settings.Slice_Order(Dummy_n)),settings.RF_Phase(Mode_n)+settings.rf_spoiling_phases(N_imaging_RF1+N_imaging_RF2),settings.IT_RF_Time,Mag_Track,settings); % Imaging RF
                    
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TR,settings.kg,settings.Diff_co,0,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation no spoiling
                    
                    % Dummy scans not stored
                    for i = 1:settings.train_spoils(IT_n)
                        P = epg_grad(P); % spoiling
                    end
                    if settings.man_spoil ==1
                        P = [0 0 P(3,1)]'; % manual spoiling
                    end
                    Cumulative_Time = Cumulative_Time + settings.IT_TR;
                end
                
                % Relaxation between last dummy train and next TR
                T_relax = settings.TR - settings.TD2 - 0.5*settings.Segment_Sizes(1).*settings.IT_TR; % Time for Relaxation between two trains(s)
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,T_relax,settings.kg,settings.Diff_co,1,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation + spoiling
                Cumulative_Time = Cumulative_Time + T_relax;
                
            end
        end % End of Dummy Scans
        
        % Start of actual scans
        for Mode_n = 1:size(settings.Tx_FA_map,3) % repeat for different modes of mTx array
            for PE2_n = 1:settings.Scan_Size(2)
                Slice_n = settings.Slice_Order(PE2_n);
                for Segment_n = 1:size(settings.Segment_Sizes,2)
                    
                    % Pre-pulse
                    N_prep_RF = N_prep_RF +1;
                    [P,~,Mag_Track] = epg_rf(P,PP_FAs(Slice_n),settings.RF_Phase(Mode_n),settings.PP_RF_Time,Mag_Track,settings);
                    Cumulative_Time = Cumulative_Time + settings.PP_RF_Time;
                    
                    for i = 1:settings.prep_spoils(N_prep_RF)
                        P = epg_grad(P); % spoiling
                    end
                    if settings.man_spoil ==1
                        P = [0 0 P(3,1)]'; % manual spoiling
                    end
                    
                    
                    
                    % Image train 1
                    for IT_n = sum(settings.Segment_Sizes(1:Segment_n - 1))+1 : sum(settings.Segment_Sizes(1:Segment_n))
                        N_imaging_RF1 = N_imaging_RF1 +1; % increment rf counter
                        [P,~,Mag_Track] = epg_rf(P,IT1_FAs(Slice_n),settings.RF_Phase(Mode_n)+settings.rf_spoiling_phases(N_imaging_RF1+N_imaging_RF2),settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                        
                        [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TE,settings.kg,settings.Diff_co,0,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation no spoiling
                        Train1(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*settings.rf_spoiling_phases(N_imaging_RF1+N_imaging_RF2)); % Store Phase-Demodulated 'signal' F+0
                        
                        [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,IT_rTE,settings.kg,settings.Diff_co,0,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation no spoiling
                        
                        for i = 1:settings.train_spoils(IT_n)
                            P = epg_grad(P); % spoiling
                        end
                        if settings.man_spoil == 1
                            P = [0 0 P(3,1)]'; % manual spoiling
                        end
                        Cumulative_Time = Cumulative_Time + settings.IT_TR;
                    end
                    
                    % Relaxation between dummy trains
                    T_relax = settings.TD2 - settings.PP_RF_Time - 1.5*settings.Segment_Sizes(1).*settings.IT_TR; % Time for Relaxation between two trains(s)
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,T_relax,settings.kg,settings.Diff_co,1,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation + spoiling
                    Cumulative_Time = Cumulative_Time + T_relax;
                    
                    % Image train 2
                    for IT_n = sum(settings.Segment_Sizes(1:Segment_n - 1))+1 : sum(settings.Segment_Sizes(1:Segment_n))
                        N_imaging_RF2 = N_imaging_RF2 +1; % increment rf counter
                        [P,~,Mag_Track] = epg_rf(P,IT2_FAs(Slice_n),settings.RF_Phase(Mode_n)+settings.rf_spoiling_phases(N_imaging_RF1+N_imaging_RF2),settings.IT_RF_Time,Mag_Track,settings); % Imaging RF
                        
                        [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TE,settings.kg,settings.Diff_co,0,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation no spoiling
                        Train2(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*settings.rf_spoiling_phases(N_imaging_RF1+N_imaging_RF2)); % Store Phase-Demodulated 'signal' F+0
                        
                        [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,IT_rTE,settings.kg,settings.Diff_co,0,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation no spoiling
                        
                        % Dummy Scans not stored
                        for i = 1:settings.train_spoils(IT_n)
                            P = epg_grad(P); % spoiling
                        end
                        if settings.man_spoil ==1
                            P = [0 0 P(3,1)]'; % manual spoiling
                        end
                        Cumulative_Time = Cumulative_Time + settings.IT_TR;
                    end
                    
                    % Relaxation between last dummy train and next TR
                    T_relax = settings.TR - settings.TD2 - 0.5*settings.Segment_Sizes(1).*settings.IT_TR; % Time for Relaxation between two trains(s)
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,T_relax,settings.kg,settings.Diff_co,1,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation + spoiling
                    Cumulative_Time = Cumulative_Time + T_relax;
                end
            end
        end
        %%% END OF SEQUENCE %%%
    end
simulation_results.IT1 = IT1;
simulation_results.IT2 = IT2;
simulation_results.N_imaging_RF1 = N_imaging_RF1;
simulation_results.N_imaging_RF2 = N_imaging_RF2;
simulation_results.N_prep_RF = N_prep_RF;
simulation_results.Cumulative_Time = Cumulative_Time;
simulation_results.Total_Energy = Total_Energy;
simulation_results.Average_10s_Power = Average_10s_Power;
end