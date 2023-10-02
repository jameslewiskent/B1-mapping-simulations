function [settings] = dream_settings(settings)
settings.title_string = 'e) DREAM';

settings.nomIT_FA = 15*pi/180; % Nominal Image Train Flip Angle in Radians
settings.IT_RF_Type = 'wSINC';
settings.IT_RF_Time = 0.6e-3;
settings.IT_RF_TBP = 3; % Asymetric windowed sinc pulse
settings.IT_RF_Pulse = Get_RF_Pulse(settings.nomIT_FA,settings.IT_RF_Type,settings.IT_RF_Time,settings.IT_RF_TBP,settings.Ref_Voltage);

settings.nomPP_FA = 55*pi/180; % Nominal Pre-pulse Flip Angle in Radians
settings.PP_RF_Type = 'wSINC';
settings.PP_RF_Time = 1.2e-3; % Time for Preparation RF pulse
settings.PP_RF_TBP = 4;
settings.PP_RF_Pulse = Get_RF_Pulse(settings.nomPP_FA,settings.PP_RF_Type,settings.PP_RF_Time,settings.PP_RF_TBP,settings.Ref_Voltage);

plot_rf_pulses(settings)

settings.prep_spoils = 10; % Number of unit gradients to move through after STEAM preparation (spoiling)
settings.noadd = 0; % = 1 to NOT add any higher-order states to spoilers. This speeds up simulations, compromises accuracy!
settings.kg = 50e3; % k-space traversal due to gradient (rad/m) for diffusion/flow
settings.rf_spoiling_phases = cumsum((0:settings.Scan_Size(1)*settings.Scan_Size(2))*50*(pi./180)); % Imaging RF phases (rad) (for RF spoiling)
settings.TR = 3.3e-3; % readout TR
settings.TE_STE = 1.2e-3;
settings.TE_FID = 1.9e-3; 
settings.Ts = settings.TE_STE + settings.TE_FID; % Ts, Time between 2 STEAM preparation pulses (s)
settings.Td = 9e-3; % Time between first preparation pulse and imaging pulse (s)

settings.Segment_Sizes = settings.Scan_Size(1); % No segmentation for DREAM

settings.Slice_Order_Type = 'OddThenEven';

if settings.UseSyntheticData == 0
    settings.RF_Phase = zeros(1,size(settings.Tx_FA_map,3)); % Otherwise it is set in the synthetic data loop
end
end

