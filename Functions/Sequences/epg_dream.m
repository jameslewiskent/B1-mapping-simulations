function [Train1,Train2,Cumulative_Time,N_imaging_RF,N_prep_RF,Mag_Track] = epg_dream(T1,IT_FAs,PP_FAs,RF_Phase,Velocity,Angle,Diff_co,settings)
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
%
% Sequence Function

P = [0 0 1]';% Z0 = 1 (Equilibrium)
Mag_Track(:,1) = [P(3,1);0];
Train1 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3));
Train2 = zeros(settings.Scan_Size(1),settings.Scan_Size(2),size(settings.Tx_FA_map,3));

Cumulative_Time = 0;
N_prep_RF = 0;
N_imaging_RF = 0;
RF_Spoil_Phase = 0; RF_Spoil_Phase_Incr = 0;

if strcmp(settings.MSor3D,'3D')
    %%% START OF 3D SEQUENCE %%%
    for Mode_n = 1:size(settings.Tx_FA_map,3) % repeat for different modes of mTx array
        
        % STEAM PREPARATION
        N_prep_RF = N_prep_RF + 1;
        [P,~,Mag_Track] = epg_rf(P,PP_FAs(Mode_n,settings.Centre_Slice),RF_Phase(Mode_n),settings.PP_RF_Time,Mag_Track,settings); % 1st STEAM RF alpha pulse
        
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
        [P,~,Mag_Track] = epg_rf(P,PP_FAs(Mode_n,settings.Centre_Slice),RF_Phase(Mode_n),settings.PP_RF_Time,Mag_Track,settings); % 2nd STEAM RF alpha pulse
        
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
                RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
                RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
                [P,~,Mag_Track] = epg_rf(P,IT_FAs(Mode_n,Slice_n),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings); % Imaging RF
                
                % STE FIRST FID SECOND SCHEME (Fig 1.a Nehkre, K. 2012)
                Cumulative_Time = Cumulative_Time + settings.TE1;
                [P,~,~,Mag_Track] =  epg_grelax(P,T1,settings.T2,settings.TE1,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation time = TE_STE
                
                P = epg_grad(P); % Adds positive "unit" Gm gradient
                
                Train1(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*RF_Spoil_Phase); % Store Phase-Demodulated 'signal' ; %P(2,2); % record STE 'signal'
                
                Cumulative_Time = Cumulative_Time + (settings.TE2 - settings.TE1);
                [P,~,~,Mag_Track] =  epg_grelax(P,T1,settings.T2,settings.TE2 - settings.TE1,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation time
                
                P = epg_mgrad(P); % Adds negative "unit" Gm gradient
                
                Train2(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*RF_Spoil_Phase); % Store Phase-Demodulated 'signal'; % record FID 'signal'
                
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
            [P,~,Mag_Track] = epg_rf(P,PP_FAs(Mode_n,Slice_n),RF_Phase(Mode_n),settings.PP_RF_Time,Mag_Track,settings); % 1st STEAM RF alpha pulse
            
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
            [P,~,Mag_Track] = epg_rf(P,PP_FAs(Mode_n,Slice_n),RF_Phase(Mode_n),settings.PP_RF_Time,Mag_Track,settings); % 2nd STEAM RF alpha pulse
            
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
                RF_Spoil_Phase = RF_Spoil_Phase + RF_Spoil_Phase_Incr;	% Increment RF Spoiling Phase
                RF_Spoil_Phase_Incr = RF_Spoil_Phase_Incr + settings.RF_Spoiling_Increment;	% Increment increment!
                [P,~,Mag_Track] = epg_rf(P,IT_FAs(Mode_n,Slice_n),RF_Phase(Mode_n)+RF_Spoil_Phase,settings.IT_RF_Time,Mag_Track,settings); % Imaging RF
                
                % STE FIRST FID SECOND SCHEME (Fig 1.a Nehkre, K. 2012)
                Cumulative_Time = Cumulative_Time + settings.TE1;
                [P,~,~,Mag_Track] =  epg_grelax(P,T1,settings.T2,settings.TE1,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation time = TE_STE
                
                P = epg_grad(P); % Adds positive "unit" Gm gradient
                
                Train1(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*RF_Spoil_Phase); % Store Phase-Demodulated 'signal' ; %P(2,2); % record STE 'signal'
                
                Cumulative_Time = Cumulative_Time + (settings.TE2 - settings.TE1);
                [P,~,~,Mag_Track] =  epg_grelax(P,T1,settings.T2,settings.TE2 - settings.TE1,settings.kg,Diff_co,0,0,Velocity,Angle,Mag_Track,settings); % relaxation time
                
                P = epg_mgrad(P); % Adds negative "unit" Gm gradient
                
                Train2(IT_n,Slice_n,Mode_n) = P(1,1).*exp(-1i*RF_Spoil_Phase); % Store Phase-Demodulated 'signal'; % record FID 'signal'
                
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
