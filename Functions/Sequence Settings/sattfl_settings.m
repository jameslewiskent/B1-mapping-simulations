function [settings] = sattfl_settings(settings)

if strcmpi(settings.MSor3D,'default')
    settings.MSor3D = '2D'; % Set the 'default' scheme
end

if strcmpi(settings.MSor3D,'3D')
    settings.title_string = 'a) 3D SatTFL';
    settings.Dummy_ITRF = 0; % Add 'dummy' RF pulses prior to image train which help achieve steady-state
    settings.TR = 10; % Long TR time conventially corresponding to 5T1 (s)
    
    settings.nom_FA = 10*(pi/180); % Nominal Image Train Flip Angle in Radians
    settings.RF_Type = 'RECT';
    settings.RF_Time = 0.1e-3;
    settings.RF_TBP = 'NA';
    settings.RF_Pulse = Get_RF_Pulse(settings.nom_FA,settings.RF_Type,settings.RF_Time,settings.RF_TBP,settings.Ref_Voltage);
    
    settings.Slice_Order_Type = 'Linear';
    
    settings.nomPP_FA = 90*(pi/180); % Nominal Pre-pulse Flip Angle in Radians
    settings.PP_RF_Type = 'RECT';
    settings.PP_RF_Time = 0.5e-3; % Time for Preparation RF pulse
    settings.PP_RF_TBP = 'NA';
    settings.PP_RF_Pulse = Get_RF_Pulse(settings.nomPP_FA,settings.PP_RF_Type,settings.PP_RF_Time,settings.PP_RF_TBP,settings.Ref_Voltage);
elseif strcmpi(settings.MSor3D,'2D')
    settings.title_string = 'a) 2D MS SatTFL';
    settings.Dummy_ITRF = 0; % Add 'dummy' RF pulses prior to image train which help achieve steady-state
    settings.TR = 10; % Long TR time conventially corresponding to 5T1 (s)
    
    settings.Slice_Order_Type = 'OddThenEven';
    settings.Min_Delay_Time = 10e-3; % Minimum delay time between reference and prepared images
    
    settings.nom_FA = 5*(pi/180); % Nominal Image Train Flip Angle in Radians
    settings.RF_Type = 'wSINC';
    settings.RF_Time = 2e-3;
    settings.RF_TBP = 4;
    settings.RF_Pulse = Get_RF_Pulse(settings.nom_FA,settings.RF_Type,settings.RF_Time,settings.RF_TBP,settings.Ref_Voltage);
    
    settings.nomPP_FA = 90*(pi/180); % Nominal Pre-pulse Flip Angle in Radians
    settings.PP_RF_Type = 'SINC';
    settings.PP_RF_Time = 5e-3; % Time for Preparation RF pulse
    settings.PP_RF_TBP = 9;
    settings.PP_RF_Pulse = Get_RF_Pulse(settings.nomPP_FA,settings.PP_RF_Type,settings.PP_RF_Time,settings.PP_RF_TBP,settings.Ref_Voltage);
end

settings.nFreqSamples = 1;
settings.IT_TE = 1.78e-3; % Image Train Echo Time (s)
settings.IT_TR = 3.9e-3; % Image Train Repitition Time (s)
settings.Segment_Factor = 1; % Segment sequence?
settings.RF_Spoiling = 1; % RF_Spoiling on (1) or off (0)
settings.RF_Spoiling_Increment = 50*(pi/180);
settings.PE1_Reordering = 'CentricOut'; % Reordering of phase encodes, 'CentricOut', 'CentricIn', 'CentricInOut', 'LinearUp', 'LinearDown'.
settings.PE2_Reordering = 'LinearUp'; % Reordering of phase encodes, 'CentricOut', 'CentricIn', 'CentricInOut', 'LinearUp', 'LinearDown'.

settings.perform_relative_mapping = 0; % Simulate relative mapping prior to absolute maps
settings.noadd = 0; % = 1 to NOT add any higher-order states to spoilers. This speeds up simulations, compromises accuracy!
settings.man_spoil = 0; % Sets transverse magnetisation to 0 (if =1) assumes 'perfect' spoiling - useful for debugging
settings.kg = 50e3; % k-space traversal due to gradient (rad/m) for diffusion/flow
settings.SR_Module = 0; % Additional 'PERFECT' Saturation recovery module, destroys all longitudinal magnetisation prior to
settings.TD_SR = 950e-3; % Saturation recovery time delay (s)
settings.magtrack_flag = 1; % Set to 1 to prevent magnetisation tracking in sequence
settings.Lookup_T1 = 0; % T1 chosen for lookup table (s)

if settings.RF_Spoiling == 0
    disp('RF spoiling is turned off.')
    settings.RF_Spoiling_Increment = 0;
end

if settings.Lookup_T1 ~= 0
warning('Calculating the flip angle map using a lookup table. This is unconventional for satTFL which typically uses the arccosine. You turn this off in the sequence settings by setting Lookup_T1 = 0.')
end

end

