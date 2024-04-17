function  [Train1,Train2,Train3,Cumulative_Time,N_imaging_RF,N_prep_RF,Mag_Track] = epg_sandwich(T1,IT_FAs,PP_FAs,RF_Phase,Velocity,Angle,Diff_co,seq_TRs,Compartment_n,settings)
% Function to simulate image trains for Sandwich B1 mapping scheme
% Controls looping structure over different flip angles, T1s etc.
%
%   INFO:
%       Adapation to Chung et al.'s B1 mapping method. Proton density image 'sandwiched' to
%       before the pre-pulse.
%
% J. Kent. 2023. Using B.Hargreaves EPG and Bloch Functions.
%
% Sequence Function

IT_rTE = settings.IT_TR - settings.IT_TE;
NTrains = 2;
if settings.T1Corr == 1
    NTrains = 3;
end

%%% START OF SEQUENCE %%%
Train1 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3));
Train2 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3));
Train3 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3));
P = [0 0 1]'; % Start magnetisation in equilibrium
Mag_Track(:,1) = [P(3,1);0];
Cumulative_Time = 0;
N_prep_RF = 0;
N_imaging_RF = 0;
RF_Spoil_Phase = 0; RF_Spoil_Phase_Incr = 0;

% Relative mapping
if settings.perform_relative_mapping == 1
    for PE2_n = 1:settings.Scan_Size(2)
        Slice_n = settings.Slice_Order(PE2_n);
        for Segment_n = 1:size(settings.Segment_Sizes,2)
            for IT_n = sum(settings.Segment_Sizes(1:Segment_n - 1))+1 : sum(settings.Segment_Sizes(1:Segment_n))
                for Mode_n = 1:size(settings.Tx_FA_map,3) % repeat for different modes of mTx array
                    N_imaging_RF = N_imaging_RF +1;
                    RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
                    RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
                    [P,~,Mag_Track] = epg_rf(P,IT_FAs(Mode_n,Slice_n),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                    
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                    
                    % We don't currently do anything with these relative data Train1(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*settings.rf_spoiling_phases(N_imaging_RF)); % Store Phase-Demodulated 'signal' F+0
                    
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

% Absolute mapping begins
for Mode_n = 1:size(settings.Tx_FA_map,3) % repeat for different modes of mTx array
    
    if settings.Dummy_Scans ~= 0
        for Dummy_n = 1:settings.Dummy_Scans
            
            % Pretend image trains
            % Dummy RF prior to image train
            for Dummy_RF = 1:settings.Dummy_ITRF
                N_imaging_RF = N_imaging_RF +1;
                RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
                RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
                [P,~,Mag_Track] = epg_rf(P,IT_FAs(Mode_n,settings.Slice_Order(jk_mod(Dummy_n,size(settings.Slice_Order,2)))),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TR,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling (FpFmZ,T1,T2,RelaxTime,kg,Diff_co,Gon,noadd)
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
                RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
                RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
                [P,~,Mag_Track] = epg_rf(P,IT_FAs(Mode_n,settings.Slice_Order(jk_mod(Dummy_n,size(settings.Slice_Order,2)))),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TR,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                
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
            [P,~,Mag_Track] = epg_rf(P,PP_FAs(Mode_n,settings.Slice_Order(jk_mod(Dummy_n,size(settings.Slice_Order,2)))),RF_Phase(Mode_n),settings.PP_RF_Time,Mag_Track,settings);
            Cumulative_Time = Cumulative_Time + settings.PP_RF_Time;
            [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2, settings.PP_Spoiler_Duration,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
            Cumulative_Time = Cumulative_Time + settings.PP_Spoiler_Duration;
            for i = 1:settings.prep_spoils(N_prep_RF)
                P = epg_grad(P); % spoiling
            end
            if settings.man_spoil == 1
                P = [0 0 P(3,1)]'; % manual spoiling
            end
            
            % Pretend image trains
            for IT_n = 1:settings.Segment_Sizes(1)
                N_imaging_RF = N_imaging_RF +1;
                RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
                RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
                [P,~,Mag_Track] = epg_rf(P,IT_FAs(Mode_n,settings.Slice_Order(jk_mod(Dummy_n,size(settings.Slice_Order,2)))),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TR,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                
                % Dummy Scans not stored
                for i = 1:settings.train_spoils(IT_n)
                    P = epg_grad(P); % spoiling
                end
                if settings.man_spoil == 1
                    P = [0 0 P(3,1)]'; % manual spoiling
                end
                Cumulative_Time = Cumulative_Time + settings.IT_TR;
            end
            
            if settings.T1Corr == 1
                % Pretend image trains
                for IT_n = 1:settings.Segment_Sizes(1)
                    N_imaging_RF = N_imaging_RF +1;
                    RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
                    RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
                    [P,~,Mag_Track] = epg_rf(P,IT_FAs(Mode_n,settings.Slice_Order(jk_mod(Dummy_n,size(settings.Slice_Order,2)))),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                    
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TR,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                    
                    % Dummy Scans not stored
                    for i = 1:settings.train_spoils(IT_n)
                        P = epg_grad(P); % spoiling
                    end
                    if settings.man_spoil == 1
                        P = [0 0 P(3,1)]'; % manual spoiling
                    end
                    Cumulative_Time = Cumulative_Time + settings.IT_TR;
                end
            end
            
            if Compartment_n == 2 % Simulation second compartment for EF
                P(3,1) = 1;
            end
            
            % Relaxation between dummy scans
            T_relax = seq_TRs(Dummy_n) - settings.PP_RF_Time - settings.PP_Spoiler_Duration - (NTrains*settings.Segment_Sizes(1).*settings.IT_TR) - settings.Dummy_ITRF*settings.IT_TR; % Time for Relaxation (s)
            [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,T_relax,settings.kg,Diff_co,1,0,Velocity,Angle,Mag_Track,settings); % relaxation + spoiling
            Cumulative_Time = Cumulative_Time + T_relax;
        end
    end % End of Dummy Scans
    
    for PE2_n = 1:settings.Scan_Size(2)
        Slice_n = settings.Slice_Order(PE2_n);
        for Segment_n = 1:size(settings.Segment_Sizes,2)
            
            % Segment of Image Train 1
            % Dummy RF prior to image train
            for Dummy_RF = 1:settings.Dummy_ITRF
                N_imaging_RF = N_imaging_RF +1;
                RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
                RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
                [P,~,Mag_Track] = epg_rf(P,IT_FAs(Mode_n,Slice_n),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TR,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling (FpFmZ,T1,T2,RelaxTime,kg,Diff_co,Gon,noadd)
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
                RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
                RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
                [P,~,Mag_Track] = epg_rf(P,IT_FAs(Mode_n,Slice_n),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                
                Train1(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*RF_Spoil_Phase); % Store Phase-Demodulated 'signal' F+0
                
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,IT_rTE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                
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
            [P,~,Mag_Track] = epg_rf(P,PP_FAs(Mode_n,Slice_n),RF_Phase(Mode_n),settings.PP_RF_Time,Mag_Track,settings);
            Cumulative_Time = Cumulative_Time + settings.PP_RF_Time;
            [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2, settings.PP_Spoiler_Duration,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
            Cumulative_Time = Cumulative_Time + settings.PP_Spoiler_Duration;
            for i = 1:settings.prep_spoils(N_prep_RF)
                P = epg_grad(P); % spoiling
            end
            if settings.man_spoil == 1
                P = [0 0 P(3,1)]'; % manual spoiling
            end
            
            
            % First segment Image Train 2
            for IT_n = sum(settings.Segment_Sizes(1:Segment_n - 1))+1 : sum(settings.Segment_Sizes(1:Segment_n))
                N_imaging_RF = N_imaging_RF +1;
                RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
                RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
                [P,~,Mag_Track] = epg_rf(P,IT_FAs(Mode_n,Slice_n),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);
                
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                
                Train2(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*RF_Spoil_Phase); % Store Phase-Demodulated 'signal' F+0
                
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,IT_rTE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                
                for i = 1:settings.train_spoils(settings.Segment_Sizes(Segment_n))
                    P = epg_grad(P); % spoiling
                end
                
                if settings.man_spoil == 1
                    P = [0 0 P(3,1)]'; % manual spoiling
                end
                Cumulative_Time = Cumulative_Time + settings.IT_TR;
            end
            
            if settings.T1Corr == 1
                % First segment Image Train 3
                for IT_n = sum(settings.Segment_Sizes(1:Segment_n - 1))+1 : sum(settings.Segment_Sizes(1:Segment_n))
                    N_imaging_RF = N_imaging_RF +1;
                    RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
                    RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
                    [P,~,Mag_Track] = epg_rf(P,IT_FAs(Mode_n,Slice_n),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);
                    
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                    
                    Train3(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*RF_Spoil_Phase); % Store Phase-Demodulated 'signal' F+0
                    
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
            
            if Compartment_n == 2 % Simulation second compartment for EF
                P(3,1) = 1;
            end
            
            % Relaxation
            T_relax = seq_TRs(settings.Dummy_Scans+Segment_n) - settings.PP_RF_Time - settings.PP_Spoiler_Duration - (NTrains*settings.Segment_Sizes(Segment_n).*settings.IT_TR) - settings.Dummy_ITRF*settings.IT_TR; % Time for Relaxation (s)
            [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,T_relax,settings.kg,Diff_co,1,0,Velocity,Angle,Mag_Track,settings); % relaxation + spoiling
            Cumulative_Time = Cumulative_Time + T_relax;
            if settings.man_spoil == 1
                P = [0 0 P(3,1)]'; % manual spoiling
            end
            
        end
    end
end



%%% END OF SEQUENCE %%%
end
