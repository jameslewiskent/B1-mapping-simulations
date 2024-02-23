function [settings] = gre_settings(settings)
settings.title_string = 'g) Relative';

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
elseif strcmp(settings.MSor3D,'2D')
settings.nom_FA = 5*(pi/180); % Nominal Image Train Flip Angle in Radians
settings.RF_Type = 'wSINC';
settings.RF_Time = 2e-3;
settings.RF_TBP = 4;
settings.RF_Pulse = Get_RF_Pulse(settings.nom_FA,settings.RF_Type,settings.RF_Time,settings.RF_TBP,settings.Ref_Voltage);
end
settings.N_TRs = 1; % N/A for GRE sequence as implemented
settings.TR = 0; % No TR long for GRE sequence as implemented
settings.IT_TE = 1.31e-3; % Image Train Echo Time (s)
settings.IT_TR = 3.14e-3; % Image Train Repitition Time (s)
settings.Dummy_RF = 100; % Add 'dummy' RF pulses prior to reference image train which help achieve steady-state
settings.Segment_Factor = 1;
settings.RF_Spoiling = 1; % RF_Spoiling on (1) or off (0)
settings.RF_Spoiling_Increment = 50*(pi/180);
settings.PE1_Reordering = 'CentricInOut'; % Reordering of phase encodes, 'CentricOut', 'CentricIn', 'CentricInOut', 'LinearUp', 'LinearDown'.
settings.PE2_Reordering = 'LinearUp'; % Reordering of phase encodes, 'CentricOut', 'CentricIn', 'CentricInOut', 'LinearUp', 'LinearDown'.
settings.Coil_Cycle = 1; % Coil-cycle transmit channels off (0) on (1)
settings.Coil_Cycle_Order = [1,4,7,2,5,8,3,6]; % NOT used if coil_cycle is off

settings.train_spoils = 1.* ones(1,ceil(settings.Scan_Size(1)/settings.Segment_Factor));
settings.noadd = 0; % = 1 to NOT add any higher-order states to spoilers. This speeds up simulations, compromises accuracy!
settings.man_spoil = 1; % Sets transverse magnetisation to 0 (if =1) assumes 'perfect' spoiling - useful for debugging
settings.kg = 50e3; % k-space traversal due to gradient (rad/m) for diffusion/flow
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
if settings.HR_TR == 1
    disp(['Heartrate variable TR option is on, TR period non-constant. Mean TR = ',num2str(settings.TR),'s. Variance = ',num2str(settings.HR_SD),'s. Min TR = ',num2str(settings.HR_minTR),'s.'])
end

if settings.RF_Spoiling == 0
    disp('RF spoiling is turned off.')
    settings.RF_Spoiling_Increment = 0;
end

if settings.Coil_Cycle == 1 
    disp('Coil cycling is turned on.');
end


end

