function [mz] = Simulate_RF_Pulse_for_plotting(RF_Pulse,RF_Pulse_Time,Dynamic_Range,Gamma,Ref_Voltage)
%

RF_Pulse = RF_Pulse*(502/Ref_Voltage);
for iB1 = 1:size(Dynamic_Range,2)
M = [0;0;1]; % Magnetisation vector at equilibirum
gmmaHzPerG = 1e-4*Gamma;  % 1e-4 is MHz/T to Hz/G
b1 = (RF_Pulse/gmmaHzPerG)*Dynamic_Range(iB1); % excitation RF field in units of Hz
gr = zeros(1,size(RF_Pulse,2)); % gradients (NOT USED)
time_points = ones(1,size(RF_Pulse,2))*(RF_Pulse_Time./size(RF_Pulse,2)); % time (in seconds)
T1 = 1.5; T2 = 25e-3;
B0_Hz = -5000:10:5000;
[~,~,mz(:,iB1)] = bloch(b1,gr,time_points,T1,T2,B0_Hz,0,0,M(1),M(2),M(3));
end

end

