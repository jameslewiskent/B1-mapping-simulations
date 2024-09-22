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
settings.Slice_Order_Type = 'Linear';
settings.TR = 10; % (s) No TR long for coil-cycled GRE sequence as implemented (only channelwise)
elseif strcmp(settings.MSor3D,'2D')
settings.nom_FA = 5*(pi/180); % Nominal Image Train Flip Angle in Radians
settings.RF_Type = 'wSINC';
settings.RF_Time = 2e-3;
settings.RF_TBP = 4;
settings.Slice_Order_Type = 'EvenThenOdd';

if settings.Matrix_Size(2) == 1 % If single-slice
    settings.TR = 4; % (s) No TR long for coil-cycled GRE sequence as implemented (only channelwise)
else % If multi-slice
    settings.TR = 10; % (s) No TR long for coil-cycled GRE sequence as implemented (only channelwise)
end

end

settings.N_TRs = 0; % N/A for GRE sequence as implemented
settings.IT_TE = 1.31e-3; % Image Train Echo Time (s)
settings.IT_TR = 3.14e-3; % Image Train Repitition Time (s)
settings.Dummy_RF = 100; % Add 'dummy' RF pulses (# per Tx channel)
settings.Segment_Factor = 1;
settings.RF_Spoiling = 1; % RF_Spoiling on (1) or off (0)
settings.RF_Spoiling_Increment = 50*(pi/180);
settings.PE1_Reordering = 'LinearUp'; % Reordering of phase encodes, 'CentricOut', 'CentricIn', 'CentricInOut', 'LinearUp', 'LinearDown'.
settings.PE2_Reordering = 'LinearUp'; % Reordering of phase encodes, 'CentricOut', 'CentricIn', 'CentricInOut', 'LinearUp', 'LinearDown'.
settings.Coil_Cycle = 'CC'; % Coil-cycle transmit channels off (single shot 'SS' or shot-wise 'SW') on ('CC')
settings.Coil_Cycle_Order = [1,4,7,2,5,8,3,6]; % NOT used if coil_cycle is off

settings.noadd = 0; % = 1 to NOT add any higher-order states to spoilers. This speeds up simulations, compromises accuracy!
settings.man_spoil = 1; % Sets transverse magnetisation to 0 (if =1) assumes 'perfect' spoiling - useful for debugging
settings.kg = 50e3; % k-space traversal due to gradient (rad/m) for diffusion/flow
settings.magtrack_flag = 1; % Set to 1 to prevent magnetisation tracking in sequence

% Notificaitons to user
if settings.HR_TR == 1
    disp(['Heartrate variable TR option is on, TR period non-constant. Mean TR = ',num2str(settings.TR),'s. Variance = ',num2str(settings.HR_SD),'s. Min TR = ',num2str(settings.HR_minTR),'s.'])
end

if settings.RF_Spoiling == 0
    disp('RF spoiling is turned off.')
    settings.RF_Spoiling_Increment = 0;
end

end

