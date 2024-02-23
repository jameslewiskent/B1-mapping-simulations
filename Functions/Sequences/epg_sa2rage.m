function [Train1,Train2,Cumulative_Time,N_imaging_RF1,N_imaging_RF2,N_prep_RF,Mag_Track] = epg_sa2rage(T1,IT1_FAs,IT2_FAs,PP_FAs,RF_Phase,Velocity,Angle,Diff_co,settings)
% Function to simulate image trains for SA2RAGE B1 mapping scheme
% Controls looping structure over different flip angles, T1s etc.
%
%   INFO:
%
% J. Kent. 2023. Using B.Hargreaves EPG and Bloch Functions.
%
% Sequence Function

P = [0 0 1]'; % Start magnetisation in equilibrium
Mag_Track(:,1) = [P(3,1);0];
IT_rTE = settings.IT_TR - settings.IT_TE;
Train1 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3));
Train2 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3));

Cumulative_Time = 0;
N_prep_RF = 0;
N_imaging_RF1 = 0;
N_imaging_RF2 = 0;
RF_Spoil_Phase = 0; RF_Spoil_Phase_Incr = 0;
%%% START OF SEQUENCE %%%

% Relative mapping
if settings.perform_relative_mapping == 1
    for PE2_n = 1:settings.Scan_Size(2)
        Slice_n = settings.Slice_Order(PE2_n);
        for Segment_n = 1:size(settings.Segment_Sizes,2)
            for IT_n = sum(settings.Segment_Sizes(1:Segment_n - 1))+1 : sum(settings.Segment_Sizes(1:Segment_n))
                for Mode_n = 1:size(settings.Tx_FA_map,3) % repeat for different modes of mTx array
                    N_imaging_RF1 = N_imaging_RF1 +1;
                    RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
                    RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
                    [P,~,Mag_Track] = epg_rf(P,IT1_FAs(Mode_n,Slice_n),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                    
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
if settings.Dummy_Scans ~= 0
    for Dummy_n = 1:settings.Dummy_Scans
        
        if size(settings.Tx_FA_map,3) > 1
            Mode_n = size(settings.Tx_FA_map,3) - (settings.Dummy_Scans - Dummy_n); % Set mode for dummy scans to last in array
        else
            Mode_n = 1; % Set mode for dummy scans
        end
        
        % Pre-pulse for dummy scans
        N_prep_RF = N_prep_RF +1;
        [P,~,Mag_Track] = epg_rf(P,PP_FAs(Mode_n,settings.Slice_Order(jk_mod(Dummy_n,length(settings.Slice_Order)))),RF_Phase(Mode_n),settings.PP_RF_Time,Mag_Track,settings);
        Cumulative_Time = Cumulative_Time + settings.PP_RF_Time;
        
        for i = 1:settings.prep_spoils(N_prep_RF)
            P = epg_grad(P); % spoiling
        end
        if settings.man_spoil ==1
            P = [0 0 P(3,1)]'; % manual spoiling
        end
        
        % Relaxation between dummy trains
        T_relax = settings.TD1 - settings.PP_RF_Time - 0.5*settings.Segment_Sizes(1).*settings.IT_TR; % Time for Relaxation between two trains(s)
        [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,T_relax,settings.kg,Diff_co,1,0,Velocity,Angle,Mag_Track,settings); % relaxation + spoiling
        Cumulative_Time = Cumulative_Time + T_relax;
        
        % Dummy image train 1
        for IT_n = 1:settings.Segment_Sizes(1)
            N_imaging_RF1 = N_imaging_RF1 +1; % increment rf counter
            RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
            RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
            [P,~,Mag_Track] = epg_rf(P,IT1_FAs(Mode_n,settings.Slice_Order(jk_mod(Dummy_n,length(settings.Slice_Order)))),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
            
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
        
        % Relaxation between dummy trains
        T_relax = settings.TD2 - settings.TD1 - settings.Segment_Sizes(1).*settings.IT_TR; % Time for Relaxation between two trains(s)
        [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,T_relax,settings.kg,Diff_co,1,0,Velocity,Angle,Mag_Track,settings); % relaxation + spoiling
        Cumulative_Time = Cumulative_Time + T_relax;
        
        % Dummy image train 2
        for IT_n = 1:settings.Segment_Sizes(1)
            N_imaging_RF2 = N_imaging_RF2 +1; % increment rf counter
            RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
            RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
            [P,~,Mag_Track] = epg_rf(P,IT2_FAs(Mode_n,settings.Slice_Order(jk_mod(Dummy_n,length(settings.Slice_Order)))),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings); % Imaging RF
            
            [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TR,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
            
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
        [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,T_relax,settings.kg,Diff_co,1,0,Velocity,Angle,Mag_Track,settings); % relaxation + spoiling
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
            [P,~,Mag_Track] = epg_rf(P,PP_FAs(Mode_n,Slice_n),RF_Phase(Mode_n),settings.PP_RF_Time,Mag_Track,settings);
            Cumulative_Time = Cumulative_Time + settings.PP_RF_Time;
            
            for i = 1:settings.prep_spoils(N_prep_RF)
                P = epg_grad(P); % spoiling
            end
            if settings.man_spoil ==1
                P = [0 0 P(3,1)]'; % manual spoiling
            end
            
            % Relaxation
            T_relax = settings.TD1 - settings.PP_RF_Time - 0.5*settings.Segment_Sizes(Segment_n).*settings.IT_TR; % Time for Relaxation between two trains(s)
            [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,T_relax,settings.kg,Diff_co,1,0,Velocity,Angle,Mag_Track,settings); % relaxation + spoiling
            Cumulative_Time = Cumulative_Time + T_relax;
            
            % Image train 1
            for IT_n = sum(settings.Segment_Sizes(1:Segment_n - 1))+1 : sum(settings.Segment_Sizes(1:Segment_n))
                N_imaging_RF1 = N_imaging_RF1 +1; % increment rf counter
                RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
                RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
                [P,~,Mag_Track] = epg_rf(P,IT1_FAs(Mode_n,Slice_n),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings);  % Imaging RF
                
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                Train1(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*RF_Spoil_Phase); % Store Phase-Demodulated 'signal' F+0
                
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,IT_rTE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                
                for i = 1:settings.train_spoils(IT_n)
                    P = epg_grad(P); % spoiling
                end
                if settings.man_spoil == 1
                    P = [0 0 P(3,1)]'; % manual spoiling
                end
                Cumulative_Time = Cumulative_Time + settings.IT_TR;
            end
            
            % Relaxation
            T_relax = settings.TD2 - settings.TD1 - settings.Segment_Sizes(Segment_n).*settings.IT_TR; % Time for Relaxation between two trains(s)
            [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,T_relax,settings.kg,Diff_co,1,0,Velocity,Angle,Mag_Track,settings); % relaxation + spoiling
            Cumulative_Time = Cumulative_Time + T_relax;
            
            % Image train 2
            for IT_n = sum(settings.Segment_Sizes(1:Segment_n - 1))+1 : sum(settings.Segment_Sizes(1:Segment_n))
                N_imaging_RF2 = N_imaging_RF2 +1; % increment rf counter
                RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
                RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
                [P,~,Mag_Track] = epg_rf(P,IT2_FAs(Mode_n,Slice_n),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings); % Imaging RF
                
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.IT_TE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                Train2(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*RF_Spoil_Phase); % Store Phase-Demodulated 'signal' F+0
                
                [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,IT_rTE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
                
                %
                for i = 1:settings.train_spoils(IT_n)
                    P = epg_grad(P); % spoiling
                end
                if settings.man_spoil ==1
                    P = [0 0 P(3,1)]'; % manual spoiling
                end
                Cumulative_Time = Cumulative_Time + settings.IT_TR;
            end
            
            % Relaxation between last dummy train and next TR
            T_relax = settings.TR - settings.TD2 - 0.5*settings.Segment_Sizes(Segment_n).*settings.IT_TR; % Time for Relaxation between two trains(s)
            [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,T_relax,settings.kg,Diff_co,1,0,Velocity,Angle,Mag_Track,settings); % relaxation + spoiling
            Cumulative_Time = Cumulative_Time + T_relax;
        end
    end
end
%%% END OF SEQUENCE %%%

end