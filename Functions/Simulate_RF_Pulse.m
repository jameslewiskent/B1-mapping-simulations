function [Flip_Angle] = Simulate_RF_Pulse(RF_Pulse,RF_Pulse_Time,B0_Hz,T1,T2,Gamma)
% Simulate the effect of frequency on the flip angle of an RF pulse

M = [0;0;1]; % Magnetisation vector at equilibirum

gmmaHzPerG = 1e-4*Gamma;  % 1e-4 is MHz/T to Hz/G

% Hz to G = b1_Hz/Gamma
b1 = RF_Pulse/gmmaHzPerG; % excitation RF field in units of Hz
gr = zeros(1,size(RF_Pulse,2)); % gradients (NOT USED)
time_points = ones(1,size(RF_Pulse,2))*(RF_Pulse_Time./size(RF_Pulse,2)); % time (in seconds)

[~,~,mz] = bloch(b1,gr,time_points,T1,T2,B0_Hz,0,0,M(1),M(2),M(3));
Flip_Angle =  acos(mz) ; % Flip angle in radians
end