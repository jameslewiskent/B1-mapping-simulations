function [mz] = Simulate_RF_Pulse_for_plotting(RF_Pulse,RF_Pulse_Time,Dynamic_Range,Gamma,Ref_Voltage)
%

RF_Pulse = RF_Pulse*(502/Ref_Voltage);
for iB1 = 1:size(Dynamic_Range,2)
    M = [0;0;1]; % Magnetisation vector at equilibirum
    gmmaHzPerG = 1e-4*Gamma;  % 1e-4 is MHz/T to Hz/G
    b1 = (RF_Pulse/gmmaHzPerG)*Dynamic_Range(iB1); % excitation RF field in units of Hz
    gr = zeros(size(RF_Pulse,2),3); % gradients (NOT USED)
    time_points = ones(size(RF_Pulse,2),1)*(RF_Pulse_Time./size(RF_Pulse,2)); % time (in seconds)
    T1 = 1.5; T2 = 25e-3;
    B0_Hz = (-5000:10:5000)';
    
    if ispc || ismac
        [~,~,mz(:,iB1)] = bloch(b1,gr,time_points,T1,T2,B0_Hz,0,0,M(1),M(2),M(3));
    else
        [~,~,mz(:,iB1)] = blochCim(b1,gr,time_points,T1,T2,B0_Hz,0,0,1); n = 1; % Mode 4 appears broken?
        while any(isnan(mz))
            % JK not sure why but this version of Bloch simulator randomly gives
            % NaNs sometimes.
            disp(['NaNs in array. Re-running for ',iptnum2ordinal(n),' time.'])
            [~,~,mz] = blochCim(b1,gr,time_points,T1,T2,B0_Hz',Slice_Positions,0,1); n = n + 1; % Mode 4 appears broken?
        end
    end
    
end
end

