function [RF_Pulse,RF_Sample_Time,RF_Pulse_Shape] = Get_RF_Pulse(nomFA,RF_Pulse_Type,RF_Pulse_Time,TBP,Ref_Voltage)
% Get RF pulse (amplitude in volts) for a nominal FA 
RF_Samples = 1e3; % Default sample size

% Reference Pulse
Ref_RF_Samples = 500;
Ref_RF_Pulse_Time = 0.5e-3; % 0.5 ms RECT pulse achieves 90 degree FA
Ref_RF_Pulse_Shape = ones(1,Ref_RF_Samples);
Ref_Amp_Integral = abs(sum(abs(Ref_RF_Pulse_Shape).*exp(1i*angle(Ref_RF_Pulse_Shape))));
Ref_RF_Sample_Time = Ref_RF_Pulse_Time/Ref_RF_Samples; 

% Get pulse shape
if strcmpi(RF_Pulse_Type,'RECT')
    RF_Pulse_Shape = ones(1,RF_Samples); 
elseif strcmpi(RF_Pulse_Type,'SINC') || strcmpi(RF_Pulse_Type,'wSINC')
    RF_Pulse_Shape = sinc(linspace(-ceil(TBP/2),floor(TBP/2),RF_Samples));
    % Window sinc
    if strcmpi(RF_Pulse_Type,'wSINC')
        RF_Pulse_Shape = RF_Pulse_Shape.*hann(size(RF_Pulse_Shape,2))';
    end
elseif strcmpi(RF_Pulse_Type,'GAUSS')
    RF_Pulse_Shape = Get_Gauss;
elseif strcmpi(RF_Pulse_Type,'HS4')
    disp('Using hard-coded HS4')
    RF_Pulse_Shape = Get_HS4;
elseif strcmpi(RF_Pulse_Type(1:2),'HS')
    RF_Pulse_Shape = Get_HSN(str2double(RF_Pulse_Type(3)),TBP); % Scaled relative to RECT
elseif strcmpi(RF_Pulse_Type,'ARB')
    RF_Pulse_Shape = Get_ARB;
else
    error('@Get_RF_Pulse: Unknown RF pulse type. Valid options are: RECT, SINC, wSINC or HSN where N = 1...8')
end
RF_Samples = size(RF_Pulse_Shape,2); % Overwrite default value with actual number of samples
RF_Sample_Time = RF_Pulse_Time/RF_Samples;

% Calculate amplitude integral
Amp_Integral = abs(sum(abs(RF_Pulse_Shape).*exp(1i*angle(RF_Pulse_Shape))));

Amp_Integral_Scale = (Ref_Amp_Integral./Amp_Integral);

Time_Scale = (Ref_RF_Sample_Time./RF_Sample_Time);

% Scale pulse voltage based on amplitude integral
RF_Pulse = (2*nomFA/pi)*RF_Pulse_Shape*Amp_Integral_Scale*Time_Scale*Ref_Voltage;

end

