function [settings] = sattfl_settings(settings)

if strcmpi(settings.MSor3D,'default')
    settings.MSor3D = '2D'; % Set the 'default' scheme
end

if strcmpi(settings.MSor3D,'3D')
    settings.title_string = 'a) 3D SatTFL';
    settings.Dummy_ITRF = 0; % Add 'dummy' RF pulses prior to image train which help achieve steady-state
    settings.TR = 10; % Long TR time conventially corresponding to 5T1 (s)
    
    settings.nomIT_FA = 10*(pi/180); % Nominal Image Train Flip Angle in Radians
    settings.IT_RF_Type = 'RECT';
    settings.IT_RF_Time = 0.1e-3;
    settings.IT_RF_TBP = 'NA';
    settings.IT_RF_Pulse = Get_RF_Pulse(settings.nomIT_FA,settings.IT_RF_Type,settings.IT_RF_Time,settings.IT_RF_TBP,settings.Ref_Voltage);
    
    settings.Slice_Shifts = zeros(1,settings.Scan_Size(2));
    settings.PP_Shifts = zeros(1,settings.Scan_Size(2));
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
    
    settings.nomIT_FA = 5*(pi/180); % Nominal Image Train Flip Angle in Radians
    settings.IT_RF_Type = 'wSINC';
    settings.IT_RF_Time = 2e-3;
    settings.IT_RF_TBP = 4;
    settings.IT_RF_Pulse = Get_RF_Pulse(settings.nomIT_FA,settings.IT_RF_Type,settings.IT_RF_Time,settings.IT_RF_TBP,settings.Ref_Voltage);
    
    settings.PP_Shifts = settings.Slice_Shifts;
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
settings.RF_Spoiling_Increment_Checks = 20;
settings.rf_spoiling_phases = cumsum((0:(4*size(settings.Tx_FA_map,3)*settings.Scan_Size(1)*settings.Scan_Size(2) + 2*settings.Dummy_ITRF*settings.Segment_Factor*size(settings.Tx_FA_map,3) + settings.RF_Spoiling_Increment_Checks))*50).*(pi./180); % Imaging phases (rad) (for RF spoiling- need enough for optional dummy scans too)
settings.rf_spoiling_phases = settings.rf_spoiling_phases(settings.RF_Spoiling_Increment_Checks:end); % Accounts for checks run on scanner which increment RF spoiling increment
settings.PE1_Reordering = 'CentricOut'; % Reordering of phase encodes, 'CentricOut', 'CentricIn', 'CentricInOut', 'LinearUp', 'LinearDown'.
settings.PE2_Reordering = 'LinearUp'; % Reordering of phase encodes, 'CentricOut', 'CentricIn', 'CentricInOut', 'LinearUp', 'LinearDown'.

settings.perform_relative_mapping = 0; % Simulate relative mapping prior to absolute maps
settings.prep_spoils = 1.* ones(1,settings.Segment_Factor*size(settings.Tx_FA_map,3)+1); % Number of unit gradients to move through after pre-pulse (spoiling- need enough for optional dummy scans too)
settings.train_spoils = 1.* ones(1,ceil(settings.Scan_Size(1)/settings.Segment_Factor));
settings.noadd = 0; % = 1 to NOT add any higher-order states to spoilers. This speeds up simulations, compromises accuracy!
settings.man_spoil = 0; % Sets transverse magnetisation to 0 (if =1) assumes 'perfect' spoiling - useful for debugging
settings.kg = 50e3; % k-space traversal due to gradient (rad/m) for diffusion/flow
settings.SR_Module = 0; % Additional 'PERFECT' Saturation recovery module, destroys all longitudinal magnetisation prior to
settings.TD_SR = 950e-3; % Saturation recovery time delay (s)
settings.magtrack_flag = 1; % Set to 1 to prevent magnetisation tracking in sequence
settings.Lookup_T1 = 0; % T1 chosen for lookup table (s)

% Calculate number of phase encodes in each segment and check
Train_Size_Tot = settings.Scan_Size(1);
settings.Segment_Sizes = zeros(1,settings.Segment_Factor);
for Seg_n = 1:settings.Segment_Factor
    settings.Segment_Sizes(Seg_n)= ceil(Train_Size_Tot./(settings.Segment_Factor+1-Seg_n));
    Train_Size_Tot = settings.Scan_Size(1) - sum(settings.Segment_Sizes,'all');
end
if sum(settings.Segment_Sizes,'all') ~= settings.Scan_Size(1)
    error('ERROR: epg_func: Number of phase encodes in segments does not equal size of train requested.')
elseif settings.Segment_Factor ~= 1
    disp(['Number of phase encodes in image train: ',num2str(settings.Scan_Size(1)),', per segment: ',num2str(settings.Segment_Sizes)]);
end

if settings.UseSyntheticData == 0
    settings.RF_Phase = zeros(1,size(settings.Tx_FA_map,3)); % Otherwise it is set in the synthetic data loop
end

if settings.RF_Spoiling == 0
    disp('RF spoiling is turned off.')
    settings.rf_spoiling_phases = zeros(size(settings.rf_spoiling_phases));
end

if settings.MTx == 1
    disp(['Multi-transmit mode mapping is active. Simulating B1 Mapping of ', num2str(size(settings.Tx_FA_map,3)),' Transmit Mode Configurations.']);
end

if settings.Lookup_T1 ~= 0
warning('Calculating the flip angle map using a lookup table. This is unconventional for satTFL which typically uses the arcosine. You turn this off in the sequence settings by setting Lookup_T1 = 0.')
end

% Need to calculate scan duration, to ensure a minimum duration between subsequent slices
Reference_Duration = settings.Scan_Size(2)*(settings.Segment_Sizes(1)*settings.IT_TR + settings.PP_RF_Time); % Minimum duration of all reference images
if Reference_Duration > settings.TR
    settings.Compulsory_Delay_Time = settings.Min_Delay_Time;
else
    settings.Compulsory_Delay_Time = (settings.TR - Reference_Duration)/(settings.Scan_Size(2)*settings.Segment_Factor);
end

end

