% Script to handle simulation and plotting of B1 mapping sequences.
% Choose a sequence to simulate and simulation parameters.
% James Kent. 2023. Using Hargreaves Bloch and EPG Functions.
clc
close all
clearvars
cd(fileparts(matlab.desktop.editor.getActiveFilename))
addpath(genpath('Functions')); load('initialise.mat');

% ---------------------------------------------------------------------- %
% ------------------- Simulation relevant parameters ------------------- %
% ---------------------------------------------------------------------- %
% Sequence specific settings can be found in 'Functions/Sequence Settings'
settings.Dynamic_Range = 0.05:0.05:3; % Simulate a large dynamic range of B1
settings.B0_Range_Hz = 0;%[0,100,500,1000]; % B0 offsets (Hz)
settings.T1s = [0.5,1,1.5,2,2.5,3]; % Array of T1 values to simulate (s)
settings.T2 = 25e-3; % T2 (s)
settings.Repeats = 1001; % Number of noise repeats to average
settings.Noise = [NaN,60]; % Simulated Noise Levels (peak SNR in Decibels) (NaN is no noise, used for generating lookup table)
settings.Scheme = 'Sandwich'; % Simulate Chosen Pulse Sequence(s) 'SatTFL', 'Sandwich', 'DREAM', 'AFI', 'SA2RAGE' or 'ALL' which uses default sequence settings
settings.MSor3D = '3D'; % 2D or 3D
settings.Velocities = 0;%[0,0.05,0.1,0.2]; % Velocity of coherent flow (m/s)
settings.Angles = 0;%[0,0,0,0]; % Angle of coherent flow (rad)
settings.Diff_coeffs = 0;%[0,3e-9]; % Diffusion co-efficient (m^2/s) (isotropic)

settings.PE1_Resolution = 1; % 4/8, 5/8, 6/8, 7/8 or 1
settings.PE1_Partial_Fourier = 1; % 4/8, 5/8, 6/8, 7/8 or 1

% Phase encode settings for 3D schemes
settings.PE2_Resolution = 1; % 4/8, 5/8, 6/8, 7/8 or 1
settings.PE2_Partial_Fourier = 1; % 4/8, 5/8, 6/8, 7/8 or 1

settings.Matrix_Size = [32 32]; % [NPE1 1] for single-slice or [NPE1 NPE2] for volumetric sequence simulation
settings.Hz_per_Volt = 8.3667; % Calibration value 8.3667 [Hz per Volt]
settings.Ref_Voltage = 60; % Applied per channel RF reference voltage [Volts per channel]

% Not implemented yet
settings.MTx = 0; % Multi-channel mapping. Only works for synthetic data
settings.Modes = 8;
settings.Enc_Scheme = 'FE'; % Encoding Method for multiTx mapping, 'FE' or 'Invert' or 'Indiv' or 'OneOFF'
settings.UseSyntheticData = 0; % Simulate on a realistic dynamic range with synthetic data (1) or (0)
settings.Whole_body_mask = 1; % (1) to simulate whole body (0) to simulate just heart region (faster)
settings.Syn_Slice = 60; % Slice of synthetic data to use

settings.Gamma = 42.57747892e6; % [Hz per T]
settings.Slice_Shift = 2e3; % Fixed slice shift (set to 0 to not fix) [Hz]

settings.Mag_Track_FAValues = [30,60,90,120,150,180]; % Track magnetisation for specific nominal flip angles
settings.Mag_Track_T1Values = 1.5; % T1 of magnetisation vector to track
settings.Mag_Track_dt = 1e-3; % Max temporal resolution of magnetisation plot
settings.EPG_trim_threshold = 0.01; % Threshold for trimming EPG states

settings.Sum_PSF = 0; % Sum PSF if = 1 or if = 0 take only centre of PSF (Recommended leave set to 0)
settings.Lookup_Size = 1e4; % Size of lookup table

%settings.LoopValues = 0.5:0.1:2; % Additional loop which can be hijacked for various purposes e.g. TR. Saves each loop in seperate .mat file.
%settings.LoopFieldName = 'TR'; % Additional loop which can be hijacked for various purposes e.g. TR. Saves each loop in seperate .mat file.
%settings.LoopValues = [32:4:64;32:4:64]'; % Additional loop which can be hijacked for various purposes e.g. TR. Saves each loop in seperate .mat file.
%settings.LoopFieldName = 'Matrix_Size'; % Additional loop which can be hijacked for various purposes e.g. TR. Saves each loop in seperate .mat file.

settings.Format = {'PE1','PE2','Tx','DR','B0','T1','Flow','Diff','Noise','Repeats'};

% Plotting settings
plot_settings.Dyn_Range_pc = 0.1; % DR defined as linear range where the difference in the mean flip angle and standard deviation is less than settings.Dyn_Range_pc (e.g. 7%) of the nominal value
plot_settings.Dynamic_Range_Axis = 1;
plot_settings.Plot_Difference = 1;

% Sequence specific inputs are contained within scheme specific epg functions.

[results,settings] = run_and_plot(settings,plot_settings,results);
