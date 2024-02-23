function [settings] = sandwich_settings(settings)
settings.title_string = 'b) Sandwich';

if strcmpi(settings.MSor3D,'default')
    settings.MSor3D = '3D'; % Set the 'default' scheme
end

if strcmp(settings.MSor3D,'3D')
settings.nom_FA = 5*pi/180; % Nominal Image Train Flip Angle in Radians
settings.RF_Type = 'RECT';
settings.RF_Time = 0.1e-3;
settings.RF_TBP = 'NA';
settings.RF_Pulse = Get_RF_Pulse(settings.nom_FA,settings.RF_Type,settings.RF_Time,settings.RF_TBP,settings.Ref_Voltage);

settings.Slice_Shifts = zeros(1,settings.Scan_Size(2));
settings.PP_Shifts = zeros(1,settings.Scan_Size(2));
settings.Slice_Order_Type = 'Linear';

settings.nomPP_FA = 90*pi/180; % Nominal Pre-pulse Flip Angle in Radians
settings.PP_RF_Type = 'HS4'; % Hyperbolic secant
settings.PP_RF_Time = 5120e-6; % Time for Preparation RF pulse (s) (5 ms)
settings.PP_Spoiler_Duration = 9000e-6; % Time for preparation pulse spoiler

%settings.PP_RF_Type = 'Rect'; % Rect
%settings.PP_RF_Time = 500e-6; % Time for Preparation RF pulse (s)

settings.PP_RF_TBP = 15; % Time bandwidth product
settings.PP_RF_Pulse = Get_RF_Pulse(settings.nomPP_FA,settings.PP_RF_Type,settings.PP_RF_Time,settings.PP_RF_TBP,settings.Ref_Voltage);

elseif strcmp(settings.MSor3D,'2D')
settings.nom_FA = 5*(pi/180); % Nominal Image Train Flip Angle in Radians
settings.RF_Type = 'wSINC';
settings.RF_Time = 2e-3;
settings.RF_TBP = 4;
settings.RF_Pulse = Get_RF_Pulse(settings.nom_FA,settings.RF_Type,settings.RF_Time,settings.RF_TBP,settings.Ref_Voltage);

settings.PP_Shifts = settings.Slice_Shifts;
settings.Slice_Order_Type = 'OddThenEven';

settings.PP_Shifts = settings.Slice_Shifts;
settings.nomPP_FA = 90*(pi/180); % Nominal Pre-pulse Flip Angle in Radians
settings.PP_RF_Type = 'SINC';
settings.PP_RF_Time = 5e-3; % Time for Preparation RF pulse
settings.PP_Spoiler_Duration = 9000e-6; % Time for preparation pulse spoiler
settings.PP_RF_TBP = 9;
settings.PP_RF_Pulse = Get_RF_Pulse(settings.nomPP_FA,settings.PP_RF_Type,settings.PP_RF_Time,settings.PP_RF_TBP,settings.Ref_Voltage);         
end

settings.T1Corr = 0; % Additional image train for T1 correction (experimental modification)
settings.IT_TE = 0.75e-3;%1.78e-3; % Image Train Echo Time (s)
settings.IT_TR = 2.5e-3;%3.9e-3; % Image Train Repitition Time (s)
settings.TR = 1; % Long TR time (s)
settings.Dummy_Scans = 2; % Add 'dummy' scans which help achieve steady-state
settings.Dummy_ITRF = 0; % Add 'dummy' RF pulses prior to reference image train which help achieve steady-state (not prior to prepared image train)
settings.Segment_Factor = 1;
settings.RF_Spoiling = 1; % RF_Spoiling on (1) or off (0)
settings.RF_Spoiling_Increment = 50*(pi/180);
settings.PE1_Reordering = 'CentricInOut'; % Reordering of phase encodes, 'CentricOut', 'CentricIn', 'CentricInOut', 'LinearUp', 'LinearDown'.
settings.PE2_Reordering = 'LinearUp'; % Reordering of phase encodes, 'CentricOut', 'CentricIn', 'CentricInOut', 'LinearUp', 'LinearDown'.

settings.perform_relative_mapping = 0; % Simulate relative mapping prior to absolute maps
settings.prep_spoils = 1.* ones(1,settings.Dummy_Scans+settings.Segment_Factor*settings.Modes*settings.Scan_Size(2)+1); %2.^(0:(Dummy_Scans+Segment_Factor*size(Tx_FA_map,3))); % Number of unit gradients to move through after pre-pulse (spoiling- need enough for optional dummy scans too)
settings.train_spoils = 1.* ones(1,ceil(settings.Scan_Size(1)/settings.Segment_Factor));
settings.noadd = 0; % = 1 to NOT add any higher-order states to spoilers. This speeds up simulations, compromises accuracy!
settings.man_spoil = 1; % Sets transverse magnetisation to 0 (if =1) assumes 'perfect' spoiling - useful for debugging
settings.kg = 50e3; % k-space traversal due to gradient (rad/m) for diffusion/flow
settings.HR_TR = 0; % Vary the TR as if it were ECG gated to a heartrate, 0 - off, 1 - on. Essentially, not strict TR but range e.g. a 1 second TR may be simulated as 1.2-1.4 s within the simulation. Vary randomly normally distributed.
settings.HR_SD = 0.1; % Standard deviation in variable heartrate TR
settings.HR_minTR = 0.7; % Minimum TR when 'gating'. Will skip to the next TR if below this limit leading to a longer TR.
settings.magtrack_flag = 1; % Set to 1 to prevent magnetisation tracking in sequence
settings.Lookup_T1 = 1.5; % T1 chosen for lookup table (s)
% Calculate number of phase encodes in each segment and check
Train_Size_Tot = settings.Scan_Size(1);
settings.Segment_Sizes = zeros(1,settings.Segment_Factor);
for Seg_n = 1:settings.Segment_Factor
    settings.Segment_Sizes(Seg_n)= ceil(Train_Size_Tot./(settings.Segment_Factor+1-Seg_n));
    Train_Size_Tot =  settings.Scan_Size(1) - sum(settings.Segment_Sizes,'all');
end
settings.N_TRs = (settings.Dummy_Scans + size(settings.Segment_Sizes,2)*settings.Scan_Size(2))*settings.Modes;

if sum(settings.Segment_Sizes,'all') ~= settings.Scan_Size(1)
    error('ERROR: epg_func: Number of phase encodes in segments does not equal size of train requested.')
elseif settings.Segment_Factor ~= 1
    disp(['Number of phase encodes in image train: ',num2str(settings.Scan_Size(1)),', per segment: ',num2str(settings.Segment_Sizes)]);
end

% Notifications to user
if settings.HR_TR == 1
    disp(['Heartrate variable TR option is on, TR period non-constant. Mean TR = ',num2str(settings.TR),'s. Variance = ',num2str(settings.HR_SD),'s. Min TR = ',num2str(settings.HR_minTR),'s.'])
end

if settings.RF_Spoiling == 0
    disp('RF spoiling is turned off.')
    settings.RF_Spoiling_Increment = 0;
end

end

