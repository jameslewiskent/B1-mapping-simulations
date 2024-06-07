function [settings] = afi_settings(settings)
settings.title_string = 'd) AFI';

if strcmpi(settings.MSor3D,'default')
    settings.MSor3D = '3D'; % Set the 'default' scheme
end

if strcmp(settings.MSor3D,'3D')
settings.nom_FA = 60*pi/180;
settings.RF_Type = 'RECT'; % Check this! Sinc Pulse TBP 3khz
settings.RF_Time = 0.5e-3;
settings.RF_TBP = 'NA';
settings.RF_Pulse = Get_RF_Pulse(settings.nom_FA,settings.RF_Type,settings.RF_Time,settings.RF_TBP,settings.Ref_Voltage);

settings.Slice_Order_Type = 'Linear';

elseif strcmp(settings.MSor3D,'2D')
error('There is no 2D Multislice AFI defined!')
end

settings.TR1 = 20e-3; %
settings.TR2 = 100e-3; %
settings.TE = 5e-3; % Image Train Echo Time (s)
settings.Dummy_Scans = 20; % Number of 'dummy' scans to achieve steady state
settings.RF_Spoiling = 1; % RF_Spoiling on (1) or off (0)
settings.RF_Spoiling_Increment = 117*(pi/180);
settings.PE1_Reordering = 'LinearUp'; % Reordering of phase encodes, 'CentricOut', 'CentricIn', 'CentricInOut', 'LinearUp', 'LinearDown'.
settings.PE2_Reordering = 'LinearUp'; % Reordering of phase encodes, 'CentricOut', 'CentricIn', 'CentricInOut', 'LinearUp', 'LinearDown'.
settings.num_spoils = 1; % Number of unit gradient spoils after each TR period
settings.man_spoil = 0; % Sets transverse magnetisation to 0 (if =1) assumes 'perfect' spoiling - useful for debugging
settings.kg = 50e3; % k-space traversal due to gradient (rad/m) for diffusion/flow
settings.Lookup_T1 = 1.5; % T1 chosen for lookup table (s)

settings.Segment_Factor = 1; % No segmentation for AFI

if settings.RF_Spoiling == 0
    disp('RF spoiling is turned off.')
    settings.RF_Spoiling_Increment = 0;
end
end

