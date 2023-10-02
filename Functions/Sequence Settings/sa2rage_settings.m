function [settings] = sa2rage_settings(settings)
settings.title_string = 'c) SA2RAGE';

settings.nomIT1_FA = 4*pi/180; % Nominal Image Train 1 Flip Angle in Radians
settings.nomIT2_FA = 11*pi/180; % Nominal Image Train 2 Flip Angle in Radians
settings.IT_RF_Type = 'RECT';
settings.IT_RF_Time = 0.1e-3;
settings.IT_RF_TBP = 'NA';
settings.IT1_RF_Pulse = Get_RF_Pulse(settings.nomIT1_FA,settings.IT_RF_Type,settings.IT_RF_Time,settings.IT_RF_TBP,settings.Ref_Voltage);
settings.IT2_RF_Pulse = Get_RF_Pulse(settings.nomIT2_FA,settings.IT_RF_Type,settings.IT_RF_Time,settings.IT_RF_TBP,settings.Ref_Voltage);

settings.nomPP_FA = 90*pi/180; % Nominal Saturation Flip Angle in Radians
settings.PP_RF_Type = 'RECT';
settings.PP_RF_Time = 0.5e-3; % Time for Preparation RF pulse
settings.PP_RF_TBP = 'NA'; % Time bandwidth product
settings.PP_RF_Pulse = Get_RF_Pulse(settings.nomPP_FA,settings.PP_RF_Type,settings.PP_RF_Time,settings.PP_RF_TBP,settings.Ref_Voltage);

settings.Slice_Order_Type = 'Linear';

settings.TD2 = 1800e-3; % Time delay 2 (s) (TD1 is unused as block immediantly follows saturation)
settings.IT_TE = 1.13e-3; % Image Train Echo Time (s)
settings.IT_TR = 2.8e-3; % Image Train Repitition Time (s)
settings.TR = 2.4; % Long TR time (s)
settings.Dummy_Scans = 4; % Add 'dummy' scans which help achieve steady-state
settings.Segment_Factor = 1;
settings.RF_Spoiling = 1; % RF_Spoiling on (1) or off (0)
settings.RF_Spoiling_Increment_Checks = 20;
settings.rf_spoiling_phases = cumsum((0:(2*settings.Dummy_Scans*settings.Scan_Size(1)*settings.Scan_Size(2)*settings.Segment_Factor*size(settings.Tx_FA_map,3) + 2*size(settings.Tx_FA_map,3)*settings.Scan_Size(1)*settings.Scan_Size(2) + 1*settings.Dummy_Scans + 1*settings.Segment_Factor*size(settings.Tx_FA_map,3) + settings.RF_Spoiling_Increment_Checks))*50).*(pi./180); % Imaging phases (rad) (for RF spoiling- need enough for optional dummy scans too)
settings.rf_spoiling_phases = settings.rf_spoiling_phases(settings.RF_Spoiling_Increment_Checks:end);
settings.PE1_Reordering = 'CentricOut'; % Reordering of phase encodes, 'CentricOut', 'CentricIn', 'CentricInOut', 'LinearUp', 'LinearDown'.
settings.PE2_Reordering = 'CentricOut'; % Reordering of phase encodes, 'CentricOut', 'CentricIn', 'CentricInOut', 'LinearUp', 'LinearDown'.

settings.prep_spoils = 1.* ones(1,settings.Dummy_Scans+settings.Segment_Factor*size(settings.Tx_FA_map,3)*settings.Scan_Size(2)+1); %2.^(0:(Dummy_Scans+Segment_Factor*size(Tx_FA_map,3))); % Number of unit gradients to move through after pre-pulse (spoiling- need enough for optional dummy scans too)
settings.train_spoils = 1.* ones(1,ceil(settings.Scan_Size(1)/settings.Segment_Factor));
settings.noadd = 0; % = 1 to NOT add any higher-order states to spoilers. This speeds up simulations, compromises accuracy!
settings.man_spoil = 0; % Sets transverse magnetisation to 0 (if =1) assumes 'perfect' spoiling - useful for debugging
settings.kg = 50e3; % k-space traversal due to gradient (rad/m) for diffusion/flow

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
    settings.RF_Phase = zeros(1,size(settings.Tx_FA_map,3)); % Otherwise it is set in the synthetic data loop
end

if settings.RF_Spoiling == 0
    disp('RF spoiling is turned off.')
    settings.rf_spoiling_phases = zeros(size(settings.rf_spoiling_phases));
end

if settings.MTx == 1
    disp(['Multi-transmit mode mapping is active. Simulating B1 Mapping of ', num2str(size(settings.Tx_FA_map,3)),' Transmit Mode Configuration.s']);
end
end

