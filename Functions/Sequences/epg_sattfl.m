function [Train1,Train2,Cumulative_Time,N_imaging_RF,N_prep_RF,Mag_Track] = epg_sattfl(T1,IT_FAs,PP_FAs,RF_Phase,Velocity,Angle,Diff_co,settings)
% Function to simulate image trains for SatTFL B1 mapping scheme.
% Controls looping structure over different flip angles, T1s etc.
%
%   INFO:
%       Replicates Chung, S. et al. MRM. 2010. 64(2):439-446. original
%       pre-pulse method of B1 mapping using a 5T1 TR to allow for full
%       relaxation of the longitudinal magnesation. Can turn on optional
%       saturation recovery module and alter the Delay Time (TD_SR).
%
% J. Kent. 2023. Using B.Hargreaves EPG and Bloch Functions.
%
% Sequence Function

P = [0 0 1]'; % Start magnetisation in equilibrium
Mag_Track(:,1) = [P(3,1);0];
IT_rTE = settings.IT_TR - settings.IT_TE - settings.IT_RF_Time;
Train1 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3));
Train2 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3));

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
                    [P,~,Mag_Track] = epg_rf(P,IT_FAs(Slice_n),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                    
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
if strcmp(settings.MSor3D,'3D')
    %%% START OF 3D SEQUENCE %%%
    for Mode_n = 1:size(settings.Tx_FA_map,3) % repeat for different modes of mTx array
        for PE2_n = 1:settings.Scan_Size(2)
            Slice_n = settings.Slice_Order(PE2_n);
            for Segment_n = 1:size(settings.Segment_Sizes,2)
                if settings.SR_Module == 1
                    P = [0 0 0]'; % Assume 'PERFECT' saturation, all longitudinal signal is nulled
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.TD_SR,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,Mag_Track,settings); % TD_SR relaxation no spoiling
                    Cumulative_Time = Cumulative_Time + settings.TD_SR;
                    settings.TR = 0; % TR set to zero, because it has relaxation during saturation recovery time delay instead
                end
                
                
                % Image Train 1
                % Dummy RF prior to image train
                for Dummy_RF = 1:settings.Dummy_ITRF
                    N_imaging_RF = N_imaging_RF +1;
                    RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
                    RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
                    [P,~,Mag_Track] = epg_rf(P,IT_FAs(Slice_n),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
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
                    [P,~,Mag_Track] = epg_rf(P,IT_FAs(Slice_n),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                    
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling (FpFmZ,T1,T2,RelaxTime,kg,Diff_co,Gon,noadd)
                    Train1(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*RF_Spoil_Phase); % Store Phase-Demodulated 'signal' F+0
                    
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,IT_rTE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                    for i = 1:settings.train_spoils(settings.Segment_Sizes(Segment_n))
                        P = epg_grad(P); % spoiling
                    end
                    
                    if settings.man_spoil ==1
                        P = [0 0 P(3,1)]'; % manual spoiling
                    end
                    Cumulative_Time = Cumulative_Time + settings.IT_TR;
                end
                
                if settings.SR_Module == 0
                    % Period of relaxation (this does nothing is SR_module is active)
                    T_relax = settings.TR-(settings.Segment_Sizes(Segment_n)*settings.IT_TR)-(settings.IT_TR*settings.Dummy_ITRF);
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,T_relax,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation + spoiling
                    Cumulative_Time = Cumulative_Time + T_relax;
                end
            end
        end
        
        for PE2_n = 1:settings.Scan_Size(2)
            Slice_n = settings.Slice_Order(PE2_n);
            for Segment_n = 1:size(settings.Segment_Sizes,2)
                if settings.SR_Module == 1
                    P(3,1) = 0 ; % Assume 'PERFECT' saturation, all longitudinal signal is nulled
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.TD_SR,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % TD_SR relaxation no spoiling
                    settings.TR = 0; % TR set to zero, because it has relaxation during saturation recovery time delay instead
                    Cumulative_Time = Cumulative_Time + settings.TD_SR;
                end
                
                
                % Pre-pulse
                N_prep_RF = N_prep_RF +1;
                [P,~,Mag_Track] = epg_rf(P,PP_FAs(Slice_n),RF_Phase(Mode_n),settings.PP_RF_Time,Mag_Track,settings);
                for i = 1:settings.prep_spoils
                    P = epg_grad(P); % spoiling
                end
                if settings.man_spoil ==1
                    P = [0 0 P(3,1)]'; % manual spoiling
                end
                Cumulative_Time = Cumulative_Time + settings.PP_RF_Time;
                
                
                % Image Train 2
                % Dummy RF prior to image train
                for Dummy_RF = 1:settings.Dummy_ITRF
                    N_imaging_RF = N_imaging_RF +1;
                    RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
                    RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
                    [P,~,Mag_Track] = epg_rf(P,IT_FAs(Slice_n),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
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
                    [P,~,Mag_Track] = epg_rf(P,IT_FAs(Slice_n),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                    
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation
                    Train2(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*RF_Spoil_Phase); % Store Phase-Demodulated 'signal' F+0
                    
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,IT_rTE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                    
                    for i = 1:settings.train_spoils(settings.Segment_Sizes(Segment_n))
                        P = epg_grad(P); % spoiling
                    end
                    
                    if settings.man_spoil ==1
                        P = [0 0 P(3,1)]'; % manual spoiling
                    end
                    Cumulative_Time = Cumulative_Time + settings.IT_TR;
                end
                
                if settings.SR_Module == 0
                    % Period of relaxation (this does nothing is SR_module is active)
                    T_relax = settings.TR - (settings.Segment_Sizes(Segment_n)*settings.IT_TR) - (settings.IT_TR*settings.Dummy_ITRF) - settings.PP_RF_Time;
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,T_relax,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation + spoiling
                    Cumulative_Time = Cumulative_Time + T_relax;
                end
            end
            
        end
    end
    %%% END OF 3D SEQUENCE %%%
elseif strcmp(settings.MSor3D,'2D')
    
    %%% START OF 2D SEQUENCE %%%
    for Mode_n = 1:size(settings.Tx_FA_map,3) % repeat for different modes of mTx array
        for PE2_n = 1:settings.Scan_Size(2)
            Slice_n = settings.Slice_Order(PE2_n);
            for Segment_n = 1:size(settings.Segment_Sizes,2)
                if settings.SR_Module == 1
                    P = [0 0 0]'; % Assume 'PERFECT' saturation, all longitudinal signal is nulled
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.TD_SR,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,Mag_Track,settings); % TD_SR relaxation no spoiling
                    Cumulative_Time = Cumulative_Time + settings.TD_SR;
                    settings.TR = 0; % TR set to zero, because it has relaxation during saturation recovery time delay instead
                end
                
                
                % Image Train 1
                % Dummy RF prior to image train
                for Dummy_RF = 1:settings.Dummy_ITRF
                    N_imaging_RF = N_imaging_RF +1;
                    RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
                    RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
                    [P,~,Mag_Track] = epg_rf(P,IT_FAs(Slice_n),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TR-settings.IT_RF_Time,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling (FpFmZ,T1,T2,RelaxTime,kg,Diff_co,Gon,noadd)
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
                    [P,~,Mag_Track] = epg_rf(P,IT_FAs(Slice_n),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                    
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling (FpFmZ,T1,T2,RelaxTime,kg,Diff_co,Gon,noadd)
                    Train1(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*RF_Spoil_Phase); % Store Phase-Demodulated 'signal' F+0
                    
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,IT_rTE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                    for i = 1:settings.train_spoils(settings.Segment_Sizes(Segment_n))
                        P = epg_grad(P); % spoiling
                    end
                    
                    if settings.man_spoil ==1
                        P = [0 0 P(3,1)]'; % manual spoiling
                    end
                    Cumulative_Time = Cumulative_Time + settings.IT_TR;
                end
                
                % Period of relaxation
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.Compulsory_Delay_Time,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation + spoiling
                Cumulative_Time = Cumulative_Time + settings.Compulsory_Delay_Time;
            end
            
        end
        
        
        for PE2_n = 1:settings.Scan_Size(2)
            Slice_n = settings.Slice_Order(PE2_n);
            for Segment_n = 1:size(settings.Segment_Sizes,2)
                if settings.SR_Module == 1
                    P(3,1) = 0 ; % Assume 'PERFECT' saturation, all longitudinal signal is nulled
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.TD_SR,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % TD_SR relaxation no spoiling
                    settings.TR = 0; % TR set to zero, because it has relaxation during saturation recovery time delay instead
                    Cumulative_Time = Cumulative_Time + settings.TD_SR;
                end
                
                
                % Pre-pulse
                N_prep_RF = N_prep_RF +1;
                [P,~,Mag_Track] = epg_rf(P,PP_FAs(Slice_n),RF_Phase(Mode_n),settings.PP_RF_Time,Mag_Track,settings);
                for i = 1:settings.prep_spoils
                    P = epg_grad(P); % spoiling
                end
                if settings.man_spoil ==1
                    P = [0 0 P(3,1)]'; % manual spoiling
                end
                Cumulative_Time = Cumulative_Time + settings.PP_RF_Time;
                
                
                % Image Train 2
                % Dummy RF prior to image train
                for Dummy_RF = 1:settings.Dummy_ITRF
                    N_imaging_RF = N_imaging_RF +1;
                    RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
                    RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
                    [P,~,Mag_Track] = epg_rf(P,IT_FAs(Slice_n),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TR-settings.IT_RF_Time,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling (FpFmZ,T1,T2,RelaxTime,kg,Diff_co,Gon,noadd)
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
                    [P,~,Mag_Track] = epg_rf(P,IT_FAs(Slice_n),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                    
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation
                    Train2(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*settings.rf_spoiling_phases(N_imaging_RF)); % Store Phase-Demodulated 'signal' F+0
                    
                    [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,IT_rTE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                    
                    for i = 1:settings.train_spoils(settings.Segment_Sizes(Segment_n))
                        P = epg_grad(P); % spoiling
                    end
                    
                    if settings.man_spoil ==1
                        P = [0 0 P(3,1)]'; % manual spoiling
                    end
                    Cumulative_Time = Cumulative_Time + settings.IT_TR;
                end
                
                %                     if settings.SR_Module == 0
                %                         % Period of relaxation (this does nothing is SR_module is active)
                %                         T_relax = settings.TR - (settings.Segment_Sizes(Segment_n)*settings.IT_TR) - (settings.IT_TR*settings.Dummy_ITRF) - settings.PP_RF_Time;
                %                         [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,T_relax,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation + spoiling
                %                         Cumulative_Time = Cumulative_Time + T_relax;
                %                     end
                % Period of relaxation
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.Compulsory_Delay_Time,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation + spoiling
                Cumulative_Time = Cumulative_Time + settings.Compulsory_Delay_Time;
            end
        end
    end
    %%% END OF 2D SEQUENCE %%%
end
end