function [settings] = afi_settings(settings)
settings.title_string = 'd) AFI';

if strcmpi(settings.MSor3D,'default')
    settings.MSor3D = '3D'; % Set the 'default' scheme
end

if strcmp(settings.MSor3D,'3D')
settings.nomFA = 60*pi/180;
settings.RF_Type = 'RECT'; % Check this! Sinc Pulse TBP 3khz
settings.RF_Time = 0.5e-3;
settings.RF_TBP = 'NA';
settings.RF_Pulse = Get_RF_Pulse(settings.nomFA,settings.RF_Type,settings.RF_Time,settings.RF_TBP,settings.Ref_Voltage);

settings.Slice_Shifts = zeros(1,settings.Scan_Size(2));
settings.PP_Shifts = zeros(1,settings.Scan_Size(2));
settings.Slice_Order_Type = 'Linear';

elseif strcmp(settings.MSor3D,'2D')
error('There is no 2D Multislice AFI defined!')
end

settings.TR1 = 20e-3; %
settings.TR2 = 100e-3; %
settings.TE = 5e-3; % Image Train Echo Time (s)
settings.Dummy_Scans = 20; % Number of 'dummy' scans to achieve steady state
settings.RF_Spoiling = 1; % RF_Spoiling on (1) or off (0)
settings.rf_spoiling_phases = cumsum((0:2*(settings.Scan_Size(1)*settings.Scan_Size(2)+settings.Dummy_Scans))*117*(pi./180)); % Imaging phases
settings.PE1_Reordering = 'LinearUp'; % Reordering of phase encodes, 'CentricOut', 'CentricIn', 'CentricInOut', 'LinearUp', 'LinearDown'.
settings.PE2_Reordering = 'LinearUp'; % Reordering of phase encodes, 'CentricOut', 'CentricIn', 'CentricInOut', 'LinearUp', 'LinearDown'.
settings.num_spoils = 1; % Number of unit gradient spoils after each TR period
settings.man_spoil = 0; % Sets transverse magnetisation to 0 (if =1) assumes 'perfect' spoiling - useful for debugging
settings.kg = 50e3; % k-space traversal due to gradient (rad/m) for diffusion/flow

settings.Segment_Sizes = settings.Scan_Size(1); % No segmentation for AFI

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
end

