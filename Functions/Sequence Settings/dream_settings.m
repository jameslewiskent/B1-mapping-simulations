function [settings] = dream_settings(settings)

if strcmpi(settings.MSor3D,'default')
    settings.MSor3D = '3D'; % Set the 'default' scheme
end

if strcmp(settings.MSor3D,'3D')
settings.title_string = 'e) 3DREAM';

settings.nom_FA = 6*pi/180; % Nominal Image Train Flip Angle in Radians
settings.RF_Type = 'RECT';
settings.RF_Time = 0.2e-3;
settings.RF_TBP = 'NA'; % Asymetric windowed sinc pulse
settings.RF_Pulse = Get_RF_Pulse(settings.nom_FA,settings.RF_Type,settings.RF_Time,settings.RF_TBP,settings.Ref_Voltage);

settings.Slice_Shifts = zeros(1,settings.Scan_Size(2));
settings.PP_Shifts = zeros(1,settings.Scan_Size(2));
settings.Slice_Order_Type = 'CentreOut';

settings.nomPP_FA = 60*pi/180; % Nominal Pre-pulse Flip Angle in Radians
settings.PP_RF_Type = 'RECT';
settings.PP_RF_Time = 0.4e-3; % Time for Preparation RF pulse
settings.PP_RF_TBP = 'NA';
settings.PP_RF_Pulse = Get_RF_Pulse(settings.nomPP_FA,settings.PP_RF_Type,settings.PP_RF_Time,settings.PP_RF_TBP,settings.Ref_Voltage);

settings.Neighbouring_Slice_Delay_Time = 0; % Compulsory delay between neighbouring slices
settings.TR = 3.06e-3; % readout TR
settings.TE1 = 1.2e-3;
settings.TE2 = 1.9e-3; 
settings.Ts = settings.TE1 + settings.TE2; % Ts, Time between 2 STEAM preparation pulses (s)
settings.Td = 8e-3; % Time between first preparation pulse and imaging pulse (s)
settings.Echo_Order = 'STEFirst'; % 'FIDFirst' or 'STEFirst'
elseif strcmp(settings.MSor3D,'2D')
settings.title_string = 'e) 2D MS DREAM';

settings.nom_FA = 10*pi/180; % Nominal Image Train Flip Angle in Radians
settings.RF_Type = 'wSINC';
settings.RF_Time = 0.7e-3;
settings.RF_TBP = 3; % Asymetric windowed sinc pulse
settings.RF_Pulse = Get_RF_Pulse(settings.nom_FA,settings.RF_Type,settings.RF_Time,settings.RF_TBP,settings.Ref_Voltage);

settings.PP_Shifts = settings.Slice_Shifts;
settings.Slice_Order_Type = 'OddThenEven';

settings.nomPP_FA = 40*pi/180; % Nominal Pre-pulse Flip Angle in Radians
settings.PP_RF_Type = 'wSINC';
settings.PP_RF_Time = 0.7e-3; % Time for Preparation RF pulse
settings.PP_RF_TBP = 4;
settings.PP_RF_Pulse = Get_RF_Pulse(settings.nomPP_FA,settings.PP_RF_Type,settings.PP_RF_Time,settings.PP_RF_TBP,settings.Ref_Voltage);

settings.Neighbouring_Slice_Delay_Time = 5; % 5 second compulsory delay between neighbouring slices
settings.TR = 3.3e-3; % readout TR
settings.TE1 = 1.2e-3;
settings.TE2 = 1.9e-3; 
settings.Ts = settings.TE1 + settings.TE2; % Ts, Time between 2 STEAM preparation pulses (s)
settings.Td = 8e-3; % Time between first preparation pulse and imaging pulse (s)
settings.Echo_Order = 'STEFirst'; % 'FIDFirst' or 'STEFirst'
end
settings.RF_Spoiling = 1; % RF_Spoiling on (1) or off (0)
settings.RF_Spoiling_Increment = 50*(pi/180);
settings.prep_spoils = 10; % Number of unit gradients to move through after STEAM preparation (spoiling)
settings.noadd = 0; % = 1 to NOT add any higher-order states to spoilers. This speeds up simulations, compromises accuracy!
settings.kg = 50e3; % k-space traversal due to gradient (rad/m) for diffusion/flow
settings.Lookup_T1 = 1.5; % T1 chosen for lookup table (s)

settings.Segment_Sizes = settings.Scan_Size(1); % No segmentation for DREAM

if settings.UseSyntheticData == 0
    settings.RF_Phase = zeros(1,size(settings.Tx_FA_map,3)); % Otherwise it is set in the synthetic data loop
end

if settings.TE2 < settings.TE1
    error('First echo time is longer than the second!')    
end

if settings.RF_Spoiling == 0
    disp('RF spoiling is turned off.')
    settings.RF_Spoiling_Increment = 0;
end

% Need to calculate time between neihgbouring slices
% Calculate slice ordering
settings = Calc_Slice_Order(settings);
Num_Slices_Between_Adjacent = find(settings.Slice_Order == 2) - find(settings.Slice_Order == 1) -1; % Find how many slices are between neighbouring slices
if isempty(Num_Slices_Between_Adjacent)
    Num_Slices_Between_Adjacent = 0;
end
Slice_Duration = (settings.Segment_Sizes(1)*settings.IT_TR + settings.Td); % Minimum duration of slice acquisition
% Delay time split across all slices
settings.Compulsory_Delay_Time = max(0,(settings.Neighbouring_Slice_Delay_Time - (Num_Slices_Between_Adjacent*Slice_Duration))/(Num_Slices_Between_Adjacent+1));
end

