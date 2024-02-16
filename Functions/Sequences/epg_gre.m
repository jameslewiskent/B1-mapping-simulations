function [simulation_results,settings] = epg_gre(settings)
% Function to simulate image trains for relative GRE mapping scheme
% Controls looping structure over different flip angles, T1s etc.
%
%   INFO:
%       Adapation to Chung et al.'s B1 mapping method. Proton density image 'sandwiched' to
%       before the pre-pulse.
%
% J. Kent. 2023. Using B.Hargreaves EPG and Bloch Functions.

% End of checks and command output, following is control of loops through T1/TR/Ratio values etc.

NRepeats = settings.Repeats; N_imaging_RF = 0;
IT1 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.IT_Tx_FA_map,3),length(settings.Dynamic_Range),length(settings.B0_Range_Hz),length(settings.T1s),length(settings.Velocities),length(settings.Diff_coeffs),1,NRepeats);
seq_TRs = zeros(settings.Dummy_Scans+size(settings.Segment_Sizes,2)*settings.Scan_Size(2),length(settings.Dynamic_Range),length(settings.B0_Range_Hz),length(settings.T1s),NRepeats);
IT_FAs = zeros(size(settings.IT_Tx_FA_map,3),settings.Scan_Size(2));
Total_Energy = zeros(1,length(settings.Dynamic_Range));
for Dynamic_Range_n = 1:size(settings.Dynamic_Range,2)
    for T1_n = 1:length(settings.T1s)
        T1 = settings.T1s(T1_n);
        for B0_n = 1:length(settings.B0_Range_Hz)
            for Mode_n = 1:size(settings.IT_Tx_FA_map,3)
                IT_FAs(Mode_n,:) = Simulate_RF_Pulse(settings.Hz_per_Volt*settings.Dynamic_Range(Mode_n,Dynamic_Range_n)*settings.IT_RF_Pulse,settings.IT_RF_Time,settings.B0_Range_Hz(B0_n)+settings.Slice_Shifts,T1,settings.T2,settings.Gamma); % Output FA in radians
            end
            settings.RF_Phase = angle(conj(settings.Dynamic_Range(:,Dynamic_Range_n)));
            if all(IT_FAs == 0)
                break
            end
            
            for Flow_n = 1:length(settings.Velocities)
                Velocity = settings.Velocities(Flow_n);
                Angle = settings.Angles(Flow_n);
                for Diff_n = 1:length(settings.Diff_coeffs)
                    Diff_co = settings.Diff_coeffs(Diff_n);
                    for Repeat_n = 1:NRepeats
                        settings = Check_Mag_Track_Flag(Dynamic_Range_n,T1_n,B0_n,Flow_n,Diff_n,settings);
                        [IT1(:,:,:,Dynamic_Range_n,B0_n,T1_n,Flow_n,Diff_n,1,Repeat_n),Cumulative_Time,N_imaging_RF,seq_TRs(:,Dynamic_Range_n,B0_n,T1_n,Repeat_n),Mag_Track] = GRE_Seq(T1,IT_FAs,Velocity,Angle,Diff_co,settings); % Simulate sequence
                        if any(settings.Mag_Track_Flags == 1)
                            simulation_results.Mag_Track{settings.Mag_Track_Flags == 1} = Mag_Track;
                        end
                    end
                end
            end
        end
    end
    if settings.verbose == 1
        disp(['Simulations ',num2str(100*Dynamic_Range_n./length(settings.Dynamic_Range),'%.2f'),'% complete.']);
    end
    Imaging_RF_Energy = Calc_RF_Energy(settings.Dynamic_Range(Dynamic_Range_n)*settings.IT_RF_Pulse,settings.IT_RF_Time); % Calculate Imaging RF pulse energy
    Total_Energy(Dynamic_Range_n) = N_imaging_RF*Imaging_RF_Energy;
end
Average_10s_Power = 10*Total_Energy./Cumulative_Time;

% Sequence Function
    function [Train1,Cumulative_Time,N_imaging_RF,seq_TRs,Mag_Track] = GRE_Seq(T1,IT_FAs,Velocity,Angle,Diff_co,settings)
        % Function runs inside relevent epg_SCHEME function and simulates a single set of image trains for specified parameters.
        
        P = [0 0 1]'; % Start magnetisation in equilibrium
        Mag_Track(:,1) = [P(3,1);0];
        IT_rTE = settings.IT_TR - settings.IT_TE;
        Train1 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3));
        
        Cumulative_Time = 0;
        N_imaging_RF = 0;
        
        seq_TRs = settings.TR.*ones(1,settings.Dummy_Scans+size(settings.Segment_Sizes,2)*settings.Scan_Size(2));
        
        
        %%% START OF SEQUENCE %%%
        % Dummy RF
        if settings.Coil_Cycle == 0
            for Mode_n = 1:size(settings.IT_Tx_FA_map,3) % repeat for different modes of mTx array
                for Dummy_RF = 1:settings.Dummy_ITRF
                    for PE2_n = 1:settings.Scan_Size(2)
                        Slice_n = settings.Slice_Order(PE2_n);
                        for IT_n = sum(settings.Segment_Sizes(1:Segment_n - 1))+1 : sum(settings.Segment_Sizes(1:Segment_n))
                            N_imaging_RF = N_imaging_RF +1;
                            [P,~,Mag_Track] = epg_rf(P,IT_FAs(Mode_n,Slice_n),settings.RF_Phase(Mode_n)+settings.rf_spoiling_phases(N_imaging_RF),settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                            [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TR,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling (FpFmZ,T1,T2,RelaxTime,kg,Diff_co,Gon,noadd)
                            for i = 1:settings.train_spoils(settings.Segment_Sizes(Segment_n))
                                P = epg_grad(P); % spoiling
                            end
                            if settings.man_spoil ==1
                                P = [0 0 P(3,1)]'; % manual spoiling
                            end
                            Cumulative_Time = Cumulative_Time + settings.IT_TR;
                        end
                    end
                end
            end
        elseif settings.Coil_Cycle == 1
            for Dummy_RF = 1:settings.Dummy_ITRF
                for PE2_n = 1:settings.Scan_Size(2)
                    Slice_n = settings.Slice_Order(PE2_n);
                    for IT_n = sum(settings.Segment_Sizes(1:Segment_n - 1))+1 : sum(settings.Segment_Sizes(1:Segment_n))
                        for Mode_n = 1:size(settings.IT_Tx_FA_map,3) % repeat for different modes of mTx array
                            N_imaging_RF = N_imaging_RF +1;
                            [P,~,Mag_Track] = epg_rf(P,IT_FAs(settings.Coil_Cycle_Order(Mode_n),Slice_n),settings.RF_Phase(Mode_n)+settings.rf_spoiling_phases(N_imaging_RF),settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                            [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TR,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling (FpFmZ,T1,T2,RelaxTime,kg,Diff_co,Gon,noadd)
                            for i = 1:settings.train_spoils(settings.Segment_Sizes(Segment_n))
                                P = epg_grad(P); % spoiling
                            end
                            if settings.man_spoil ==1
                                P = [0 0 P(3,1)]'; % manual spoiling
                            end
                            Cumulative_Time = Cumulative_Time + settings.IT_TR;
                        end
                    end
                end
            end
        end
        
        % Relative mapping
        if settings.Coil_Cycle == 0
            for Mode_n = 1:size(settings.IT_Tx_FA_map,3) % repeat for different modes of mTx array
                for PE2_n = 1:settings.Scan_Size(2)
                    Slice_n = settings.Slice_Order(PE2_n);
                    for Segment_n = 1:size(settings.Segment_Sizes,2)
                        % Segment of Image Train
                        for IT_n = sum(settings.Segment_Sizes(1:Segment_n - 1))+1 : sum(settings.Segment_Sizes(1:Segment_n))
                            N_imaging_RF = N_imaging_RF +1;
                            [P,~,Mag_Track] = epg_rf(P,IT_FAs(settings.Coil_Cycle_Order(Mode_n),Slice_n),settings.RF_Phase(Mode_n)+settings.rf_spoiling_phases(N_imaging_RF),settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                            
                            [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                            
                            Train1(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*settings.rf_spoiling_phases(N_imaging_RF)); % Store Phase-Demodulated 'signal' F+0
                            
                            [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,IT_rTE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                            
                            for i = 1:settings.train_spoils(settings.Segment_Sizes(Segment_n))
                                P = epg_grad(P); % spoiling
                            end
                            
                            if settings.man_spoil == 1
                                P = [0 0 P(3,1)]'; % manual spoiling
                            end
                            Cumulative_Time = Cumulative_Time + settings.IT_TR;
                        end
                    end
                end
            end
        elseif settings.Coil_Cycle == 1
            for PE2_n = 1:settings.Scan_Size(2)
                Slice_n = settings.Slice_Order(PE2_n);
                for Segment_n = 1:size(settings.Segment_Sizes,2)
                    % Segment of Image Train
                    for IT_n = sum(settings.Segment_Sizes(1:Segment_n - 1))+1 : sum(settings.Segment_Sizes(1:Segment_n))
                        for Mode_n = 1:size(settings.IT_Tx_FA_map,3) % repeat for different modes of mTx array
                            
                            N_imaging_RF = N_imaging_RF +1;
                            [P,~,Mag_Track] = epg_rf(P,IT_FAs(settings.Coil_Cycle_Order(Mode_n),Slice_n),settings.RF_Phase(Mode_n)+settings.rf_spoiling_phases(N_imaging_RF),settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                            
                            [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                            
                            Train1(IT_n,Slice_n,settings.Coil_Cycle_Order(Mode_n)) = P(1,1).*exp(-1i*settings.rf_spoiling_phases(N_imaging_RF)); % Store Phase-Demodulated 'signal' F+0
                            
                            [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,IT_rTE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                            
                            for i = 1:settings.train_spoils(settings.Segment_Sizes(Segment_n))
                                P = epg_grad(P); % spoiling
                            end
                            
                            if settings.man_spoil == 1
                                P = [0 0 P(3,1)]'; % manual spoiling
                            end
                            Cumulative_Time = Cumulative_Time + settings.IT_TR;
                        end
                    end
                end
            end
        end
        
        %%% END OF SEQUENCE %%%
    end
simulation_results.IT1 = IT1;
simulation_results.N_imaging_RF = N_imaging_RF;
simulation_results.Cumulative_Time = Cumulative_Time;
simulation_results.Total_Energy = Total_Energy;
simulation_results.seq_TRs = seq_TRs;
simulation_results.Average_10s_Power = Average_10s_Power;
end
