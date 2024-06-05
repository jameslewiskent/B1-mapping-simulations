% Script to handle simulation and plotting of B1 mapping sequences.
% Choose a sequence to simulate and simulation parameters.
% Sequence specific settings can be found in 'Functions/Sequence Settings'
% James Kent. 2024. Using Hargreaves Bloch and EPG Functions.
clc
close all
clearvars
cd(fileparts(matlab.desktop.editor.getActiveFilename))
addpath(genpath('Functions')); load('initialise.mat');

% ---------------------------------------------------------------------- %
% ----------------------- Simulation Parameters ------------------------ %
% ---------------------------------------------------------------------- %
settings.MSor3D = '3D'; % 2D or 3D or default
settings.Scheme = 'GRE'; % Simulate Chosen Pulse Sequence(s) 'SatTFL', 'Sandwich', 'DREAM', 'AFI', 'SA2RAGE' or 'ALL' which uses default sequence settings
settings.Dynamic_Range = 0:0.01:3; % Simulate a large dynamic range of B1
settings.B0_Range_Hz = 0; %[0,100,500,1000]; % B0 offsets (Hz)
settings.T1s = 0.5:0.1:4;%[0.5,1,1.5,2,2.5,3]; % Array of T1 values to simulate (s)
settings.T2 = 25e-3; % T2 (s)
settings.Repeats = 1; % Number of noise repeats to average
settings.Noise = [NaN,60]; % Simulated Noise Levels (peak SNR in Decibels) (NaN is no noise, used for generating lookup table)
settings.Velocities = 0;%[0,0.05,0.1,0.2]; % Velocity of coherent flow (m/s)
settings.Angles = 0;%[0,0,0,0]; % Angle of coherent flow (rad)
settings.Diff_coeffs = 0;%[0,3e-9]; % Diffusion co-efficient (m^2/s) (isotropic)
settings.Matrix_Size = [32 1]; % [NPE1 1] for single-slice or [NPE1 NPE2] for volumetric sequence simulation
settings.PE1_Resolution = 1; % 4/8, 5/8, 6/8, 7/8 or 1
settings.PE1_Partial_Fourier = 1; % 4/8, 5/8, 6/8, 7/8 or 1
settings.PE2_Resolution = 1; % 4/8, 5/8, 6/8, 7/8 or 1
settings.PE2_Partial_Fourier = 1; % 4/8, 5/8, 6/8, 7/8 or 1

% ---------------------------------------------------------------------- %
% ---------------- Synthetic Body Simulation Parameters ---------------- %
% ---------------------------------------------------------------------- %
settings.UseSyntheticData = 1; % Simulate on a realistic dynamic range with synthetic data (1) or (0). (Replaces Dynamic_Range and T1s with synthetic body values)
settings.Modes = 8; % Multi-channel modes. (1) is CP. Only works for synthetic data with defined Tx sensitivites
settings.Enc_Scheme = 'Indiv'; % Encoding scheme for mTx mapping, 'FE' or 'Invert' or 'Indiv' or 'OneOFF' or 'OneXpc' where X = 0 to 100
settings.T1_blood = 2.29; % T1 of arterial blood [s]
settings.T1_heart = 1.925; % T1 of heart [s]
settings.Global_T1 = 1; % Use a fixed T1 over the entire body (given by a single value in settings.T1s)
settings.Whole_body_mask = 1; % (1) to simulate whole body (0) to simulate just heart region (1) which is quicker
settings.Syn_Slice = 60; % Slice of synthetic data to use
settings.Mag_Track_SynInd = [92,40;34,25;11,59;8,103;16,143;50,157;95,139;104,75;82,89;62,92]; % Track magnetisation for voxel in heart [76,100] 
settings.Ejection_Fraction = 0; % Ejection fraction - mixes two compartments magnetisation each TR with specified ratio (only works if implemented within epg_seq.m)

% ---------------------------------------------------------------------- %
% ------------------------- Advanced Settings -------------------------- %
% ---------------------------------------------------------------------- %
settings.Gamma = 42.57747892e6; % [Hz per T]
settings.Slice_Thickness = 1; % Slice Thickness [cm]
settings.Hz_per_Volt = 8.3667; % Calibration value 8.3667 [Hz per Volt]
settings.Ref_Voltage = 60; % Applied per channel RF reference voltage [Volts per channel]
settings.Mag_Track_FAValues = [30,90,150]; % Track magnetisation for specific nominal flip angles
settings.Mag_Track_T1Values = 1.5; % T1 of magnetisation vector to track
settings.Mag_Track_dt = 1e-3; % Max temporal resolution of magnetisation plot
settings.EPG_trim_threshold = 0.01; % Threshold for trimming EPG states
settings.Calc_FWHM = 0; % Calculate the FWHM of the PSF
settings.Sum_PSF = 0; % Sum PSF if = 1 or if = 0 take only centre of PSF (Recommended leave set to 0)
settings.Lookup_Size = 1e4; % Size of lookup table
settings.Use_Previous_Lookup = 0; % Use previous lookup table (1) or generate new (default) (0)
settings.verbose = 1;

% Additional loop which can be hijacked for various purposes e.g. TR. Saves each loop in seperate .mat file.
% Uncomment required loop field name and values
% settings.LoopFieldName = 'HR_TR'; settings.LoopValues = [0,1]';
% settings.LoopFieldName = 'Ejection_Fraction'; settings.LoopValues = [0:0.01:1]';
% settings.LoopFieldName = 'TR'; settings.LoopValues = [0.5:0.05:2]'; 
settings.LoopFieldName = 'Coil_Cycle'; settings.LoopValues = [0,1]';
%settings.LoopFieldName2 = 'Coil_Cycle_Order'; settings.LoopValues2 = [1:8;1,4,7,2,5,8,3,6];
%settings.LoopFieldName2 = 'nom_FA'; settings.LoopValues2 = ([1:15]*pi/180)';

settings.Format = {'PE1','PE2','Tx','DR','B0','T1','Flow','Diff','Noise','Repeats'};

% ---------------------------------------------------------------------- %
% ---------------------------- Plot Settings --------------------------- %
% ---------------------------------------------------------------------- %
plot_settings.Dyn_Range_pc = 0.1; % DR defined as linear range where the difference in the mean flip angle and standard deviation is less than settings.Dyn_Range_pc (e.g. 7%) of the nominal value
plot_settings.Dynamic_Range_Axis = 0;
plot_settings.Plot_Difference = 0;
plot_settings.Show_Dyn_Range = 1; % Plots green box over dynamic range

tic
[results,settings] = run_and_plot(settings,plot_settings,results);
toc