function [Flip_Angle] = Simulate_RF_Pulse(RF_Pulse,RF_Pulse_Time,B0_Hz,Slice_Positions,T1,T2,Gamma)
% Simulate the effect of frequency on the flip angle of an RF pulse

M = [0;0;1]; % Magnetisation vector at equilibirum

gmmaHzPerG = 1e-4*Gamma;  % 1e-4 is MHz/T to Hz/G

% Hz to G = b1_Hz/Gamma
b1 = (RF_Pulse/gmmaHzPerG)'; % excitation RF field in units of Hz
time_points = ones(size(RF_Pulse,2),1)*(RF_Pulse_Time./size(RF_Pulse,2)); % time (in seconds)

if any(Slice_Positions(:,1) ~= 0)
    gr = [ones(size(RF_Pulse,2),1).*0.5,zeros(size(RF_Pulse,2),2)]; % gradients 1,2 or 3-dimensional gradient in G/cm. 5mT/m = 0.5G/cm    
else
    gr = zeros(size(RF_Pulse,2),3); % gradients (NOT USED) 1,2 or 3-dimensional gradient in G/cm.
end

if ispc || ismac
    [~,~,mz] = bloch(b1,gr,time_points,T1,T2,B0_Hz',Slice_Positions,0,M(1),M(2),M(3));
else
    [~,~,mz] = blochCim(b1,gr,time_points,T1,T2,B0_Hz',Slice_Positions,0,1); n = 1; % Mode 4 appears broken?
    while any(isnan(mz))
        % JK not sure why but this version of Bloch simulator randomly gives
        % NaNs sometimes.
        disp(['NaNs in array. Re-running for ',iptnum2ordinal(n),' time.'])
        [~,~,mz] = blochCim(b1,gr,time_points,T1,T2,B0_Hz',Slice_Positions,0,1); n = n + 1; % Mode 4 appears broken?
    end
end

Flip_Angle =  acos(mz); % Flip angle in radians
end