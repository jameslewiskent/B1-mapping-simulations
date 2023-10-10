function [RF_Pulse,RF_Sample_Time] = Get_RF_Pulse(nomFA,RF_Pulse_Type,RF_Pulse_Time,TBP,Ref_Voltage)
% Get RF pulse (amplitude in volts) for a nominal 90 deg FA 
RF_Samples = 1e3;

if nomFA == 0
    error('Nominal flip angle cannot be 0.')
end

% Get pulse shape
if strcmpi(RF_Pulse_Type,'RECT')
    RF_Pulse_Shape = ones(1,RF_Samples); 
elseif strcmp(RF_Pulse_Type,'SINC') || strcmp(RF_Pulse_Type,'wSINC')
    RF_Pulse_Shape = sinc(linspace(-ceil(TBP/2),floor(TBP/2),RF_Samples));
    % Window sinc
    if strcmp(RF_Pulse_Type,'wSINC')
        RF_Pulse_Shape = RF_Pulse_Shape.*hann(size(RF_Pulse_Shape,2))';
    end
elseif strcmpi(RF_Pulse_Type,'GAUSS')
    error('GAUSS RF PULSE NOT YET DEFINED');
elseif strcmpi(RF_Pulse_Type(1:2),'HS')
    HSN = str2double(RF_Pulse_Type(3));
    RF_Pulse_Shape = Get_HSN(HSN,TBP); % Scaled relative to RECT
    RF_Samples = 1024; % Hard set in Get_HSN
else
    error('@Get_RF_Pulse: Unknown RF pulse type. Valid options are: RECT, SINC, wSINC or HSN where N = 1...8')
end

% Scale pulse
RF_Pulse = (2*nomFA/pi)*RF_Pulse_Shape*Scale_RF_Pulse(RF_Pulse_Shape,RF_Pulse_Time,Ref_Voltage);

RF_Sample_Time = RF_Pulse_Time/RF_Samples;
end

