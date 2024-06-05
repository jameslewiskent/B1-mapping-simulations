function [Train1,Cumulative_Time,N_imaging_RF,Mag_Track]  = epg_gre(T1,IT_FAs,RF_Phase,Velocity,Angle,Diff_co,settings)
% Function to simulate image trains for relative GRE mapping scheme
% Controls looping structure over different flip angles, T1s etc.
%
% J. Kent. 2023.

% Sequence Function

P = [0 0 1]'; % Start magnetisation in equilibrium
Mag_Track(:,1) = [P(3,1);0];
IT_rTE = settings.IT_TR - settings.IT_TE;
Train1 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3));

Cumulative_Time = 0;
N_imaging_RF = 0;
RF_Spoil_Phase = 0; RF_Spoil_Phase_Incr = 0;

%%% START OF SEQUENCE %%%

% Dummy RF
if settings.Coil_Cycle == 0
    for Mode_n = 1:size(settings.Tx_FA_map,3) % repeat for different modes of mTx array
        if N_imaging_RF == settings.Dummy_RF
            break
        end
        for PE2_n = 1:settings.Scan_Size(2)
            Slice_n = settings.Slice_Order(PE2_n);
            if N_imaging_RF == settings.Dummy_RF
                break
            end
            for Segment_n = 1:size(settings.Segment_Sizes,2)
                if N_imaging_RF == settings.Dummy_RF
                    break
                end
                for IT_n = sum(settings.Segment_Sizes(1:Segment_n - 1))+1 : sum(settings.Segment_Sizes(1:Segment_n))
                    N_imaging_RF = N_imaging_RF +1;
                    RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
                    RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
                    [P,~,Mag_Track] = epg_rf(P,IT_FAs(settings.Coil_Cycle_Order(Mode_n),Slice_n),RF_Phase(settings.Coil_Cycle_Order(Mode_n))+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TR,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling (FpFmZ,T1,T2,RelaxTime,kg,Diff_co,Gon,noadd)
                    for i = 1:settings.train_spoils(settings.Segment_Sizes(Segment_n))
                        P = epg_grad(P); % spoiling
                    end
                    if settings.man_spoil ==1
                        P = [0 0 P(3,1)]'; % manual spoiling
                    end
                    Cumulative_Time = Cumulative_Time + settings.IT_TR;
                    
                    if N_imaging_RF == settings.Dummy_RF
                        break
                    end
                end
            end
            
            % Relaxation (if any) between channels
            if settings.TR >= (settings.Scan_Size(1)*settings.IT_TR)
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.TR - (settings.Scan_Size(1)*settings.IT_TR),settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling (FpFmZ,T1,T2,RelaxTime,kg,Diff_co,Gon,noadd)
            end
            
        end
    end
elseif settings.Coil_Cycle == 1
    for PE2_n = 1:settings.Scan_Size(2)
        Slice_n = settings.Slice_Order(PE2_n);
        if N_imaging_RF == settings.Dummy_RF
            break
        end
        for Segment_n = 1:size(settings.Segment_Sizes,2)
            if N_imaging_RF == settings.Dummy_RF
                break
            end
            for IT_n = sum(settings.Segment_Sizes(1:Segment_n - 1))+1 : sum(settings.Segment_Sizes(1:Segment_n))
                if N_imaging_RF == settings.Dummy_RF
                    break
                end
                for Mode_n = 1:size(settings.Tx_FA_map,3) % repeat for different modes of mTx array
                    N_imaging_RF = N_imaging_RF + 1;
                    RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
                    RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
                    [P,~,Mag_Track] = epg_rf(P,IT_FAs(settings.Coil_Cycle_Order(Mode_n),Slice_n),RF_Phase(settings.Coil_Cycle_Order(Mode_n))+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TR,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling (FpFmZ,T1,T2,RelaxTime,kg,Diff_co,Gon,noadd)
                    for i = 1:settings.train_spoils(settings.Segment_Sizes(Segment_n))
                        P = epg_grad(P); % spoiling
                    end
                    if settings.man_spoil ==1
                        P = [0 0 P(3,1)]'; % manual spoiling
                    end
                    Cumulative_Time = Cumulative_Time + settings.IT_TR;
                    
                    if N_imaging_RF == settings.Dummy_RF
                        break
                    end
                end
            end
        end
    end
end

% Relative mapping
if settings.Coil_Cycle == 0
    for Mode_n = 1:size(settings.Tx_FA_map,3) % repeat for different modes of mTx array
        for PE2_n = 1:settings.Scan_Size(2)
            Slice_n = settings.Slice_Order(PE2_n);
            for Segment_n = 1:size(settings.Segment_Sizes,2)
                % Segment of Image Train
                for IT_n = sum(settings.Segment_Sizes(1:Segment_n - 1))+1 : sum(settings.Segment_Sizes(1:Segment_n))
                    N_imaging_RF = N_imaging_RF +1;
                    RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
                    RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
                    [P,~,Mag_Track] = epg_rf(P,IT_FAs(settings.Coil_Cycle_Order(Mode_n),Slice_n),RF_Phase(settings.Coil_Cycle_Order(Mode_n))+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                    
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                    
                    Train1(IT_n,Slice_n,settings.Coil_Cycle_Order(Mode_n)) = P(1,1).*exp(-1i*RF_Spoil_Phase); % Store Phase-Demodulated 'signal' F+0
                    
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
            
            % Relaxation (if any) between channels
            if settings.TR >= (settings.Scan_Size(1)*settings.IT_TR)
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.TR - (settings.Scan_Size(1)*settings.IT_TR),settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling (FpFmZ,T1,T2,RelaxTime,kg,Diff_co,Gon,noadd)
            end
            
        end
    end
elseif settings.Coil_Cycle == 1
    for PE2_n = 1:settings.Scan_Size(2)
        Slice_n = settings.Slice_Order(PE2_n);
        for Segment_n = 1:size(settings.Segment_Sizes,2)
            % Segment of Image Train
            for IT_n = sum(settings.Segment_Sizes(1:Segment_n - 1))+1 : sum(settings.Segment_Sizes(1:Segment_n))
                for Mode_n = 1:size(settings.Tx_FA_map,3) % repeat for different modes of mTx array
                    N_imaging_RF = N_imaging_RF +1;
                    RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
                    RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
                    [P,~,Mag_Track] = epg_rf(P,IT_FAs(settings.Coil_Cycle_Order(Mode_n),Slice_n),RF_Phase(settings.Coil_Cycle_Order(Mode_n))+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                    
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                    
                    Train1(IT_n,Slice_n,settings.Coil_Cycle_Order(Mode_n)) = P(1,1).*exp(-1i*RF_Spoil_Phase); % Store Phase-Demodulated 'signal' F+0
                    
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
