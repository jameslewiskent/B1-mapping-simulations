function [settings] = sa2rage_settings(settings)
settings.title_string = 'c) SA2RAGE';

if strcmpi(settings.MSor3D,'default')
    settings.MSor3D = '3D'; % Set the 'default' scheme
end

settings.nom_FA = 4*pi/180; % Nominal Image Train 1 Flip Angle in Radians
settings.nom_FA2 = 11*pi/180; % Nominal Image Train 2 Flip Angle in Radians
settings.RF_Type = 'RECT';
settings.RF_Time = 0.1e-3;
settings.RF_TBP = 'NA';
settings.RF_Pulse = Get_RF_Pulse(settings.nom_FA,settings.RF_Type,settings.RF_Time,settings.RF_TBP,settings.Ref_Voltage);
settings.RF_Pulse2 = Get_RF_Pulse(settings.nom_FA2,settings.RF_Type,settings.RF_Time,settings.RF_TBP,settings.Ref_Voltage);

settings.Slice_Shifts = zeros(1,settings.Scan_Size(2));
settings.PP_Shifts = zeros(1,settings.Scan_Size(2));
settings.Slice_Order_Type = 'Linear';

settings.nomPP_FA = 90*pi/180; % Nominal Saturation Flip Angle in Radians
settings.PP_RF_Type = 'RECT';
settings.PP_RF_Time = 0.5e-3; % Time for Preparation RF pulse
settings.PP_RF_TBP = 'NA'; % Time bandwidth product
settings.PP_RF_Pulse = Get_RF_Pulse(settings.nomPP_FA,settings.PP_RF_Type,settings.PP_RF_Time,settings.PP_RF_TBP,settings.Ref_Voltage);

settings.TD1 = 53e-3; % Time delay 1 (s)
settings.TD2 = 1800e-3; % Time delay 2 (s) 
settings.IT_TE = 1.13e-3; % Image Train Echo Time (s)
settings.IT_TR = 2.8e-3; % Image Train Repitition Time (s)
settings.TR = 2.4; % Long TR time (s)
settings.Dummy_Scans = 4; % Add 'dummy' scans which help achieve steady-state
settings.Segment_Factor = 1;
settings.RF_Spoiling = 1; % RF_Spoiling on (1) or off (0)
settings.RF_Spoiling_Increment = 50*(pi/180);
settings.PE1_Reordering = 'CentricOut'; % Reordering of phase encodes, 'CentricOut', 'CentricIn', 'CentricInOut', 'LinearUp', 'LinearDown'.
settings.PE2_Reordering = 'CentricOut'; % Reordering of phase encodes, 'CentricOut', 'CentricIn', 'CentricInOut', 'LinearUp', 'LinearDown'.

settings.perform_relative_mapping = 0; % Simulate relative mapping prior to absolute maps
settings.prep_spoils = 1.* ones(1,settings.Dummy_Scans+settings.Segment_Factor*size(settings.Tx_FA_map,3)*settings.Scan_Size(2)+1); %2.^(0:(Dummy_Scans+Segment_Factor*size(Tx_FA_map,3))); % Number of unit gradients to move through after pre-pulse (spoiling- need enough for optional dummy scans too)
settings.train_spoils = 1.* ones(1,ceil(settings.Scan_Size(1)/settings.Segment_Factor));
settings.noadd = 0; % = 1 to NOT add any higher-order states to spoilers. This speeds up simulations, compromises accuracy!
settings.man_spoil = 0; % Sets transverse magnetisation to 0 (if =1) assumes 'perfect' spoiling - useful for debugging
settings.kg = 50e3; % k-space traversal due to gradient (rad/m) for diffusion/flow
settings.Lookup_T1 = 1.5; % T1 chosen for lookup table (s)

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
if settings.RF_Spoiling == 0
    disp('RF spoiling is turned off.')
    settings.RF_Spoiling_Increment = 0;
end

if settings.TD1 < settings.PP_RF_Time + 0.5*settings.Segment_Sizes(1)*settings.IT_TR
    settings.TD1 = settings.PP_RF_Time + 0.5*settings.Segment_Sizes(1)*settings.IT_TR;
    warning(['TD1 too low. Changing to minimum delay time = ',num2str(settings.TD1),' s.'])
end

end

