function [settings] = sandwich_settings(settings)
settings.title_string = 'b) Sandwich';

settings.nomIT_FA = 4*pi/180; % Nominal Image Train Flip Angle in Radians
settings.IT_RF_Type = 'RECT';
settings.IT_RF_Time = 0.1e-3;
settings.IT_RF_TBP = 'NA';
settings.IT_RF_Pulse = Get_RF_Pulse(settings.nomIT_FA,settings.IT_RF_Type,settings.IT_RF_Time,settings.IT_RF_TBP,settings.Ref_Voltage);

settings.PP_Shifts = zeros(1,settings.Scan_Size(2));
settings.nomPP_FA = 80*pi/180; % Nominal Pre-pulse Flip Angle in Radians
settings.PP_RF_Type = 'HS8'; % Hyperbolic secant 8
settings.PP_RF_Time = 5e-3; % Time for Preparation RF pulse (s) (5 ms)
settings.PP_RF_TBP = 15; % Time bandwidth product
settings.PP_RF_Pulse = Get_RF_Pulse(settings.nomPP_FA,settings.PP_RF_Type,settings.PP_RF_Time,settings.PP_RF_TBP,settings.Ref_Voltage);

settings.Slice_Order_Type = 'Linear';

settings.IT_TE = 1.78e-3; % Image Train Echo Time (s)
settings.IT_TR = 3.9e-3; % Image Train Repitition Time (s)
settings.TR = 1; % Long TR time (s)
settings.Dummy_Scans = 2; % Add 'dummy' scans which help achieve steady-state
settings.Dummy_ITRF = 0; % Add 'dummy' RF pulses prior to reference image train which help achieve steady-state (not prior to prepared image train)
settings.Segment_Factor = 1;
settings.RF_Spoiling = 1; % RF_Spoiling on (1) or off (0)
settings.RF_Spoiling_Increment_Checks = 20;
settings.rf_spoiling_phases = cumsum((0:(2*settings.Dummy_Scans*settings.Scan_Size(1)*settings.Scan_Size(2)*settings.Segment_Factor*size(settings.Tx_FA_map,3) + 2*size(settings.Tx_FA_map,3)*settings.Scan_Size(1)*settings.Scan_Size(2) + 1*settings.Dummy_ITRF*settings.Dummy_Scans + 1*settings.Dummy_ITRF*settings.Segment_Factor*size(settings.Tx_FA_map,3) + settings.RF_Spoiling_Increment_Checks))*50).*(pi./180); % Imaging phases (rad) (for RF spoiling- need enough for optional dummy scans too)
settings.rf_spoiling_phases = settings.rf_spoiling_phases(settings.RF_Spoiling_Increment_Checks:end); % Accounts for checks run on scanner which increment RF spoiling increment
settings.PE1_Reordering = 'CentricInOut'; % Reordering of phase encodes, 'CentricOut', 'CentricIn', 'CentricInOut', 'LinearUp', 'LinearDown'.
settings.PE2_Reordering = 'LinearUp'; % Reordering of phase encodes, 'CentricOut', 'CentricIn', 'CentricInOut', 'LinearUp', 'LinearDown'.

settings.prep_spoils = 1.* ones(1,settings.Dummy_Scans+settings.Segment_Factor*size(settings.Tx_FA_map,3)*settings.Scan_Size(2)+1); %2.^(0:(Dummy_Scans+Segment_Factor*size(Tx_FA_map,3))); % Number of unit gradients to move through after pre-pulse (spoiling- need enough for optional dummy scans too)
settings.train_spoils = 1.* ones(1,ceil(settings.Scan_Size(1)/settings.Segment_Factor));
settings.noadd = 0; % = 1 to NOT add any higher-order states to spoilers. This speeds up simulations, compromises accuracy!
settings.man_spoil = 0; % Sets transverse magnetisation to 0 (if =1) assumes 'perfect' spoiling - useful for debugging
settings.kg = 50e3; % k-space traversal due to gradient (rad/m) for diffusion/flow
settings.HR_TR = 0; % Vary the TR as if it were ECG gated to a heartrate, 0 - off, 1 - on. Essentially, not strict TR but range e.g. a 1 second TR may be simulated as 1.2-1.4 s within the simulation. Vary randomly normally distributed.
settings.HR_SD = 0.1; % Standard deviation in variable heartrate TR
settings.HR_minTR = 0.7; % Minimum TR when 'gating'. Will skip to the next TR if below this limit leading to a longer TR.
settings.Ejection_Frac = 0; % Ejection fraction on (1) or off (0) - resets a proportion of the magnetisation
settings.magtrack_flag = 1; % Set to 1 to prevent magnetisation tracking in sequence


% Calculate number of phase encodes in each segment and check
Train_Size_Tot = settings.Scan_Size(1);
settings.Segment_Sizes = zeros(1,settings.Segment_Factor);
for Seg_n = 1:settings.Segment_Factor
    settings.Segment_Sizes(Seg_n)= ceil(Train_Size_Tot./(settings.Segment_Factor+1-Seg_n));
    Train_Size_Tot =  settings.Scan_Size(1) - sum(settings.Segment_Sizes,'all');
end
if sum(settings.Segment_Sizes,'all') ~= settings.Scan_Size(1)
    error('ERROR: epg_func: Number of phase encodes in segments does not equal size of train requested.')
elseif settings.Segment_Factor ~= 1
    disp(['Number of phase encodes in image train: ',num2str(settings.Scan_Size(1)),', per segment: ',num2str(settings.Segment_Sizes)]);
end

% Notificaitons to user
if settings.UseSyntheticData == 0
    settings.Ejection_Frac = 0; % Do not use ejection fraction for non-synthetic data sims
    settings.RF_Phase = zeros(1,size(settings.Tx_FA_map,3)); % Otherwise it is set in the synethic data loop
end


if settings.HR_TR == 1
    disp(['Heartrate variable TR option is on, TR period non-constant. Mean TR = ',num2str(settings.TR),'s. Variance = ',num2str(settings.HR_SD),'s. Min TR = ',num2str(settings.HR_minTR),'s.'])
end

if settings.RF_Spoiling == 0
    disp('RF spoiling is turned off.')
    settings.rf_spoiling_phases = zeros(size(settings.rf_spoiling_phases));
end

if settings.MTx == 1
    disp(['Multi-transmit mode mapping is active. Simulating B1 Mapping of ', num2str(size(settings.Tx_FA_map,3)),' Transmit Mode Configuration.s']);
end

end

