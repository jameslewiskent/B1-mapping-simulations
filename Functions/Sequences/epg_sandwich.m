function [simulation_results,settings] = epg_sandwich(settings)
% Function to simulate image trains for Sandwich B1 mapping scheme
% Controls looping structure over different flip angles, T1s etc.
%
%   INFO:
%       Adapation to Chung et al.'s B1 mapping method. Proton density image 'sandwiched' to
%       before the pre-pulse.
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
    %                 [IT1(:,indi,indj,:),IT2(:,indi,indj,:),Cumulative_Time,Total_Energy(Dynamic_Range_n),N_imaging_RF,N_prep_RF] = Sandwich_Seq(PP_FA,IT_FA,TR,T1,T2,IT_TR,IT_TE,settings.Train_Size,rf_spoiling_phases,kg,settings.Diff_co,settings.Velocity,settings.Angle,prep_spoils,train_spoils,noadd,man_spoil,HR_TR,Dummy_Scans,Ejection_Frac,magtrack_flag,Segment_Sizes,nFreqSamples,Tx_FA_map,RF_Phase);
    %             end
    %         end
    %
    %         disp(['Simulating Image Trains. Percent Complete: ',num2str(100*indi/size(slice_T1,1),'%.2f'),'%'])
    %     end
    
else
    IT1 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3),length(settings.Dynamic_Range),length(settings.B0_Range_Hz),length(settings.T1s));
    IT2 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3),length(settings.Dynamic_Range),length(settings.B0_Range_Hz),length(settings.T1s));
    seq_TRs = zeros(settings.Dummy_Scans+size(settings.Segment_Sizes,2),length(settings.Dynamic_Range),length(settings.B0_Range_Hz),length(settings.T1s));
    Total_Energy = zeros(1,length(settings.Dynamic_Range));
    for Dynamic_Range_n = 1:length(settings.Dynamic_Range)
        for T1_n = 1:length(settings.T1s)
            for B0_n = 1:length(settings.B0_Range_Hz)
                T1 = settings.T1s(T1_n);
                IT_FAs = Simulate_RF_Pulse(settings.Hz_per_Volt*settings.Dynamic_Range(Dynamic_Range_n)*settings.IT_RF_Pulse,settings.IT_RF_Time,settings.B0_Range_Hz(B0_n)+settings.Slice_Shifts,T1,settings.T2,settings.Gamma); % Output FA in radians
                PP_FAs = Simulate_RF_Pulse(settings.Hz_per_Volt*settings.Dynamic_Range(Dynamic_Range_n)*settings.PP_RF_Pulse,settings.PP_RF_Time,settings.B0_Range_Hz(B0_n)+settings.PP_Shifts,T1,settings.T2,settings.Gamma); % Output FA in radians
                settings = Check_Mag_Track_Flag(Dynamic_Range_n,T1_n,B0_n,settings);
                [IT1(:,:,:,Dynamic_Range_n,B0_n,T1_n),IT2(:,:,:,Dynamic_Range_n,B0_n,T1_n),Cumulative_Time,N_imaging_RF,N_prep_RF,seq_TRs(:,Dynamic_Range_n,B0_n,T1_n),Mag_Track] = Sandwich_Seq(T1,IT_FAs,PP_FAs,settings); % Simulate sequence
                if any(settings.Mag_Track_Flags == 1)
                    simulation_results.Mag_Track{settings.Mag_Track_Flags == 1} = Mag_Track;
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
    function [Train1,Train2,Cumulative_Time,N_imaging_RF,N_prep_RF,seq_TRs,Mag_Track] = Sandwich_Seq(T1,IT_FAs,PP_FAs,settings)
        % Function runs inside relevent epg_SCHEME function and simulates a single set of image trains for specified parameters.
        
        P = [0 0 1]'; % Start magnetisation in equilibrium
        Mag_Track(:,1) = [P(3,1);0];
        IT_rTE = settings.IT_TR - settings.IT_TE;
        Train1 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3));
        Train2 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3));
        
        Cumulative_Time = 0;
        N_prep_RF = 0;
        N_imaging_RF = 0;
        
        if settings.HR_TR == 1 % Variable TRs simulating gating to a heart rate
            HR_TR = zeros(1,settings.Dummy_Scans+size(settings.Segment_Sizes,2));
            for TR_n = 1:settings.Dummy_Scans+size(settings.Segment_Sizes,2)
                HR_TR(TR_n) = settings.TR + settings.HR_SD*settings.TR*randn; % Variable TRs
                if HR_TR(TR_n) < settings.HR_minTR % Enforce a minimum TR by skipping a TR if too short
                    HR_TR(TR_n) = HR_TR(TR_n) + settings.TR + settings.HR_SD*settings.TR*randn;
                end
            end
            seq_TRs = HR_TR;
        else
            seq_TRs = settings.TR.*ones(1,settings.Dummy_Scans+size(settings.Segment_Sizes,2));
        end
        
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
                
                if settings.Ejection_Frac == 1
                    if T1 == 2.29 % T1 of blood (hence currently simulating blood voxel)
                        P = [0 0 1]'; % manual spoiling due to ejected/refreshed blood
                    end
                end
                
                % Pretend image trains
                % Dummy RF prior to image train
                for Dummy_RF = 1:settings.Dummy_ITRF
                    N_imaging_RF = N_imaging_RF +1;
                    [P,~,Mag_Track] = epg_rf(P,IT_FAs(settings.Slice_Order(Dummy_n)),settings.RF_Phase(Mode_n)+settings.rf_spoiling_phases(N_imaging_RF),settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TR,settings.kg,settings.Diff_co,0,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation no spoiling (FpFmZ,T1,T2,RelaxTime,kg,Diff_co,Gon,noadd)
                    for i = 1:settings.train_spoils(settings.Segment_Sizes(1))
                        P = epg_grad(P); % spoiling
                    end
                    if settings.man_spoil ==1
                        P = [0 0 P(3,1)]'; % manual spoiling
                    end
                    Cumulative_Time = Cumulative_Time + settings.IT_TR;
                end
                for IT_n = 1:settings.Segment_Sizes(1)
                    N_imaging_RF = N_imaging_RF +1;
                    
                    [P,~,Mag_Track] = epg_rf(P,IT_FAs(settings.Slice_Order(Dummy_n)),settings.RF_Phase(Mode_n)+settings.rf_spoiling_phases(N_imaging_RF),settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                    
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TR,settings.kg,settings.Diff_co,0,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation no spoiling
                    
                    % Dummy Scans not stored
                    for i = 1:settings.train_spoils(IT_n)
                        P = epg_grad(P); % spoiling
                    end
                    if settings.man_spoil ==1
                        P = [0 0 P(3,1)]'; % manual spoiling
                    end
                    
                    Cumulative_Time = Cumulative_Time + settings.IT_TR;
                end
                
                % Pre-pulse for dummy scans
                N_prep_RF = N_prep_RF +1;
                [P,~,Mag_Track] = epg_rf(P,PP_FAs(settings.Slice_Order(Dummy_n)),settings.RF_Phase(Mode_n),settings.PP_RF_Time,Mag_Track,settings);
                Cumulative_Time = Cumulative_Time + settings.PP_RF_Time;
                for i = 1:settings.prep_spoils(N_prep_RF)
                    P = epg_grad(P); % spoiling
                end
                if settings.man_spoil == 1
                    P = [0 0 P(3,1)]'; % manual spoiling
                end
                
                % Pretend image trains
                for IT_n = 1:settings.Segment_Sizes(1)
                    N_imaging_RF = N_imaging_RF +1;
                    [P,~,Mag_Track] = epg_rf(P,IT_FAs(settings.Slice_Order(Dummy_n)),settings.RF_Phase(Mode_n)+settings.rf_spoiling_phases(N_imaging_RF),settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                    
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
                
                % Relaxation between dummy scans
                T_relax = seq_TRs(1,Dummy_n) - settings.PP_RF_Time - (2*settings.Segment_Sizes(1).*settings.IT_TR) - settings.Dummy_ITRF*settings.IT_TR; % Time for Relaxation (s)
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,T_relax,settings.kg,settings.Diff_co,1,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation + spoiling
                Cumulative_Time = Cumulative_Time + T_relax;
            end
        end % End of Dummy Scans
        
        
        for Mode_n = 1:size(settings.Tx_FA_map,3) % repeat for different modes of mTx array
            if settings.Ejection_Frac == 1
                if T1 == 2.29 % T1 of blood (hence currently simulating blood voxel)
                    P = [0 0 1]'; % manual spoiling due to ejected/refreshed blood
                end
            end
            for PE2_n = 1:settings.Scan_Size(2)
                Slice_n = settings.Slice_Order(PE2_n);
                for Segment_n = 1:size(settings.Segment_Sizes,2)
                    
                    % Segment of Image Train 1
                    % Dummy RF prior to image train
                    for Dummy_RF = 1:settings.Dummy_ITRF
                        N_imaging_RF = N_imaging_RF +1;
                        [P,~,Mag_Track] = epg_rf(P,IT_FAs(Slice_n),settings.RF_Phase(Mode_n)+settings.rf_spoiling_phases(N_imaging_RF),settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                        [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TR,settings.kg,settings.Diff_co,0,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation no spoiling (FpFmZ,T1,T2,RelaxTime,kg,Diff_co,Gon,noadd)
                        for i = 1:settings.train_spoils(settings.Segment_Sizes(Segment_n))
                            P = epg_grad(P); % spoiling
                        end
                        if settings.man_spoil ==1
                            P = [0 0 P(3,1)]'; % manual spoiling
                        end
                        Cumulative_Time = Cumulative_Time + settings.IT_TR;
                    end
                    for IT_n = sum(settings.Segment_Sizes(1:Segment_n - 1))+1 : sum(settings.Segment_Sizes(1:Segment_n))
                        N_imaging_RF = N_imaging_RF +1;
                        [P,~,Mag_Track] = epg_rf(P,IT_FAs(Slice_n),settings.RF_Phase(Mode_n)+settings.rf_spoiling_phases(N_imaging_RF),settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                        
                        [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TE,settings.kg,settings.Diff_co,0,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation no spoiling
                        
                        Train1(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*settings.rf_spoiling_phases(N_imaging_RF)); % Store Phase-Demodulated 'signal' F+0
                        
                        [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,IT_rTE,settings.kg,settings.Diff_co,0,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation no spoiling
                        
                        for i = 1:settings.train_spoils(settings.Segment_Sizes(Segment_n))
                            P = epg_grad(P); % spoiling
                        end
                        
                        if settings.man_spoil == 1
                            P = [0 0 P(3,1)]'; % manual spoiling
                        end
                        Cumulative_Time = Cumulative_Time + settings.IT_TR;
                    end
                    
                    
                    % Pre-pulse
                    N_prep_RF = N_prep_RF +1;
                    [P,~,Mag_Track] = epg_rf(P,PP_FAs(Slice_n),settings.RF_Phase(Mode_n),settings.PP_RF_Time,Mag_Track,settings);
                    Cumulative_Time = Cumulative_Time + settings.PP_RF_Time;
                    for i = 1:settings.prep_spoils(N_prep_RF)
                        P = epg_grad(P); % spoiling
                    end
                    if settings.man_spoil == 1
                        P = [0 0 P(3,1)]'; % manual spoiling
                    end
                    
                    
                    % First segment Image Train 2
                    for IT_n = sum(settings.Segment_Sizes(1:Segment_n - 1))+1 : sum(settings.Segment_Sizes(1:Segment_n))
                        N_imaging_RF = N_imaging_RF +1;
                        [P,~,Mag_Track] = epg_rf(P,IT_FAs(Slice_n),settings.RF_Phase(Mode_n)+settings.rf_spoiling_phases(N_imaging_RF),settings.IT_RF_Time,Mag_Track,settings);
                        
                        [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TE,settings.kg,settings.Diff_co,0,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation no spoiling
                        
                        Train2(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*settings.rf_spoiling_phases(N_imaging_RF)); % Store Phase-Demodulated 'signal' F+0
                        
                        [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,IT_rTE,settings.kg,settings.Diff_co,0,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation no spoiling
                        
                        for i = 1:settings.train_spoils(settings.Segment_Sizes(Segment_n))
                            P = epg_grad(P); % spoiling
                        end
                        
                        if settings.man_spoil == 1
                            P = [0 0 P(3,1)]'; % manual spoiling
                        end
                        Cumulative_Time = Cumulative_Time + settings.IT_TR;
                    end
                    
                    if settings.Ejection_Frac == 1
                        if T1 == 2.29 % T1 of blood (hence currently simulating blood voxel)
                            P = [0 0 1]'; % manual spoiling due to ejected/refreshed blood
                        end
                    end
                    
                    % Relaxation
                    T_relax = seq_TRs(1,settings.Dummy_Scans+Segment_n) - settings.PP_RF_Time - (2*settings.Segment_Sizes(Segment_n).*settings.IT_TR) - settings.Dummy_ITRF*settings.IT_TR; % Time for Relaxation (s)
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,T_relax,settings.kg,settings.Diff_co,1,0,settings.Velocity,settings.Angle,Mag_Track,settings); % relaxation + spoiling
                    Cumulative_Time = Cumulative_Time + T_relax;
                    if settings.man_spoil == 1
                        P = [0 0 P(3,1)]'; % manual spoiling
                    end
                    
                end
            end
        end
        %%% END OF SEQUENCE %%%
    end
simulation_results.IT1 = IT1;
simulation_results.IT2 = IT2;
simulation_results.N_imaging_RF = N_imaging_RF;
simulation_results.N_prep_RF = N_prep_RF;
simulation_results.Cumulative_Time = Cumulative_Time;
simulation_results.Total_Energy = Total_Energy;
simulation_results.seq_TRs = seq_TRs;
simulation_results.Average_10s_Power = Average_10s_Power;
end
