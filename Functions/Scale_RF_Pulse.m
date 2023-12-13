function [Scale_Factor] = Scale_RF_Pulse(RF_Pulse_Shape,RF_Pulse_Time,Ref_Voltage)
% Scale RF pulse to achieve 90 degree flip angle relatove to rect

% M = [0;0;1]; % Magnetisation vector at equilibirum
% T1 = 2; T2 = 25e-3;
% 
% gmmaHzPerG = 1e2*42.57747892;  % 1e2 is MHz/T to Hz/G
% % Hz to G = b1_Hz/Gamma
% 
% B1_Range_Hz = 0:1:5000;
% mz = zeros(length(B1_Range_Hz),1);
% for B1_n = 1:length(B1_Range_Hz)
%     b1 = RF_Pulse_Shape*B1_Range_Hz(B1_n)/gmmaHzPerG; % excitation RF field in units of Hz
%     gr = zeros(1,size(RF_Pulse_Shape,2)); % gradients (NOT USED)
%     time_points = ones(1,size(RF_Pulse_Shape,2))*(RF_Pulse_Time./size(RF_Pulse_Shape,2)); % time (in seconds)
%     
%     [~,~,mz(B1_n,:)] = bloch(b1,gr,time_points,T1,T2,0,0,0,M(1),M(2),M(3));
% end
% indices = find(abs(mz) <= 0.01); % Finds list of indices which are close to zero
% 
% if isempty(indices) % Check an index has been found
%     warning('Unable to scale RF pulse. Retrying covering a larger B1 range.')
%     B1_Range_Hz = 0:1:50000;
%     mz = zeros(length(B1_Range_Hz),1);
%     for B1_n = 1:length(B1_Range_Hz)
%         b1 = RF_Pulse_Shape*B1_Range_Hz(B1_n)/gmmaHzPerG; % excitation RF field in units of Hz
%         gr = zeros(1,size(RF_Pulse_Shape,2)); % gradients (NOT USED)
%         time_points = ones(1,size(RF_Pulse_Shape,2))*(RF_Pulse_Time./size(RF_Pulse_Shape,2)); % time (in seconds)
%         [~,~,mz(B1_n,:)] = bloch(b1,gr,time_points,T1,T2,0,0,0,M(1),M(2),M(3));
%     end
%     indices = find(abs(mz) <= 0.01); % Finds list of indices which are close to zero
%     if isempty(indices)
%         error('Warning! Unable to scale RF pulse. Please check you are requesting an achievable pulse.')
%     end
% end
% [~,index] = min(abs(mz(indices(1:find([diff(indices)>1;1],1,'first'))))); % Find which index is first and closest to zero



Scale_Factor = B1_Range_Hz(indices(index))/(502/Ref_Voltage); % 0.5 ms RECT achieves 90 degree flip angle at 8.37 Hz at 60V
end

