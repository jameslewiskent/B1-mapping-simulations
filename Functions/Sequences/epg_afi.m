function [Train1,Train2,Cumulative_Time,N_imaging_RF,Mag_Track] = epg_afi(T1,FAs,RF_Phase,Velocity,Angle,Diff_co,settings)
% Function to simulate image trains for AFI B1 mapping scheme.
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
% J. Kent. 2023.
%
% Sequence Function

P = [0 0 1]'; % Z0 = 1 (Equilibrium)
Mag_Track(:,1) = [P(3,1);0];
rTE1 = settings.TR1 - settings.TE;
rTE2 = settings.TR2 - settings.TE;
Train1 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3));
Train2 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3));

Cumulative_Time = 0;
N_imaging_RF = 0;
RF_Spoil_Phase = 0; RF_Spoil_Phase_Incr = 0;

%%% START OF SEQUENCE %%%
for Mode_n = 1:size(settings.Tx_FA_map,3) % repeat for different modes of mTx array
    
    % Dummy repetitions to produce steady state
    for Dummy_n = 1:settings.Dummy_Scans
        N_imaging_RF = N_imaging_RF + 1;
        RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
        RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
        [P,~,Mag_Track] = epg_rf(P,FAs(Mode_n,settings.Slice_Order(jk_mod(Dummy_n,settings.Scan_Size(2)))),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.RF_Time,Mag_Track,settings);
        [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.TR1,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
        Cumulative_Time = Cumulative_Time + settings.TR1;
        for i = 1:settings.num_spoils
            P = epg_grad(P); % spoiling
        end
        if settings.man_spoil == 1
            P = [0 0 P(3,1)]'; % manual spoiling
        end
        
        N_imaging_RF = N_imaging_RF + 1;
        RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
        RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
        [P,~,Mag_Track] = epg_rf(P,FAs(Mode_n,settings.Slice_Order(jk_mod(Dummy_n,settings.Scan_Size(2)))),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.RF_Time,Mag_Track,settings);
        [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.TR2,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
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
            RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
            RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
            [P,~,Mag_Track] = epg_rf(P,FAs(Mode_n,Slice_n),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.RF_Time,Mag_Track,settings);
            
            [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.TE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
            
            Train1(TR_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*RF_Spoil_Phase); % Store Phase-Demodulated 'signal' F+0
            [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,rTE1,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
            Cumulative_Time = Cumulative_Time + settings.TR1;
            for i = 1:settings.num_spoils
                P = epg_grad(P); % spoiling
            end
            if settings.man_spoil ==1
                P = [0 0 P(3,1)]'; % manual spoiling
            end
            N_imaging_RF = N_imaging_RF +1;
            RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
            RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
            [P,~,Mag_Track] = epg_rf(P,FAs(Mode_n,Slice_n),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.RF_Time,Mag_Track,settings);
            
            [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,settings.TE,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
            
            Train2(TR_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*RF_Spoil_Phase); % Store Phase-Demodulated 'signal' F+0
            [P,~,~,Mag_Track] = epg_grelax(P,T1,settings.T2,rTE2,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation no spoiling
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

%%% END OF SEQUENCE %%%
end
