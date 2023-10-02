% Script to handle simulation and plotting of B1 mapping sequences.
% Choose a sequence to simulate and simulation parameters.
%
% James Kent. 2023. Using Hargreaves Bloch and EPG Functions.
clc
close all
clear all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
addpath(genpath('Functions')); load('initialise.mat');

% ---------------------------------------------------------------------- %
% ------------------- Simulation relevant parameters ------------------- %
% ---------------------------------------------------------------------- %
% Sequence specific settings can be found in 'Functions/Sequence Settings'

settings.Dynamic_Range = 0.05:0.05:3; % Simulate a large dynamic range of B1
settings.B0_Range_Hz = 0:500:1500; % B0 Frequency offset (Hz)
settings.T1s = [0.5,1,1.5,2,2.5,3]; % Array of T1 values to simulate (s)
settings.T2 = 25e-3; % T2 (s)
settings.Repeats = 18; % Number of noise repeats to average
settings.Noise = 60; % Simulated Noise Levels (peak SNR in Decibels)
settings.Scheme = 'AFI'; % Simulate Chosen Pulse Sequence(s) 'SatTFL','Sandwich', 'DREAM', 'AFI', 'SA2RAGE' or 'ALL'
settings.MSor3D = '3D'; % Multi-slice (MS) or 3D
settings.Velocity = 0; % Velocity of coherent flow (m/s)
settings.Angle = 0; % Angle of coherent flow (rad)
settings.Diff_co = 0; % Diffusion co-efficient (m^2/s) (isotropic)

settings.PE1_Resolution = 1; % 4/8, 5/8, 6/8, 7/8 or 1
settings.PE1_Partial_Fourier = 1; % 4/8, 5/8, 6/8, 7/8 or 1

% Phase encode settings for 3D schemes
settings.PE2_Resolution = 1; % 4/8, 5/8, 6/8, 7/8 or 1
settings.PE2_Partial_Fourier = 1; % 4/8, 5/8, 6/8, 7/8 or 1

settings.Matrix_Size = [32 32]; % [NPE1 1] for 2D single-slice or [NPE1 NPE2] for 2D multi-slice or 3D sequence simulation
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
settings.Mag_Track_T1Values = settings.T1s(3); % T1 of magnetisation vector to track
settings.Mag_Track_dt = 1e-3; % Max temporal resolution of magnetisation plot
settings.EPG_trim_threshold = 0.01; % Threshold for trimming EPG states

settings.Sum_PSF = 0; % Sum PSF if = 1 or if = 0 take only centre of PSF (Recommended leave set to 0)
settings.Lookup_T1 = 1.5; % T1 chosen for lookup table (s)

% Plotting settings
plot_settings.Dyn_Range_pc = 0.1; % DR defined as linear range where the difference in the mean flip angle and standard deviation is less than settings.Dyn_Range_pc (e.g. 7%) of the nominal value
plot_settings.Dynamic_Range_Axis = 1;
plot_settings.Plot_Difference = 1;

% Sequence specific inputs are contained within scheme specific epg
% functions.
% ---------------------------------------------------------------------- %
% -------------------------- END of user input ------------------------- %
% ---------------------------------------------------------------------- %

[results,settings] = Masterscript(settings,plot_settings);

function [results,settings] = Masterscript(settings,plot_settings)
% This is the core function which controls plotting of the data. Hands off
% to other functions to simulate data if necessary- e.g. if not already simulated with identical parameters.

% ---------------------------------------------------------------------- %
% ----------- Following code handles simulation and plotting ----------- %
% ---------------------------------------------------------------------- %
if settings.UseSyntheticData == 0
    
    Schemes = {'SatTFL','Sandwich','SA2RAGE','AFI','DREAM'};
    if strcmp(settings.Scheme,'ALL')
        settings = repmat(settings,1,length(Schemes));
        for Scheme_n = 1:length(Schemes)
            settings(Scheme_n).Scheme = Schemes{Scheme_n};
            [results(1,Scheme_n),settings(1,Scheme_n)] = run_sequence_simulations(settings(1,Scheme_n)); % Simulate and process data
        end
        
        figure('color','w'); tiledlayout('flow');
        for Scheme_n = 1:length(Schemes)
        nexttile; plot_T1_fig(results(1,Scheme_n),settings(1,Scheme_n),plot_settings); % T1 Plot
        end
        
        plot_SAR_fig(results,settings); % SAR Plot

    else
    [results,settings] = run_sequence_simulations(settings); % Simulate and process data
    figure('color','w'); plot_T1_fig(results,settings,plot_settings); % T1 Plot
    figure('color','w'); plot_B0_fig(results,settings,plot_settings); % B0 Plot
    figure('color','w'); plot_mz(results,settings); % magnetisation plot
    end
    
else
    
    % ---------------------------------------------------------------------- %
    % -----------             Simulating Synthetic Data          ----------- %
    % ---------------------------------------------------------------------- %
    
    % Read in and generate slice data if not already existing
    SyntheticDuke = load('Data\SyntheticDuke.mat'); % Reads in if current folder is Masterscript
    % Assign T1, T2 values to synthetic data, mask out non-heart regions (temp - to speed up dev)
    
    slice = SyntheticDuke.sigma(:,:,settings.Syn_Slice); % For now only consider one slice of conductivity data to mask
    Gamma = 42.57747892e6; % Gyromagnetic ratio of Hydrogen (Hz/T)
    
    % Define relaxation constants for blood & heart
    settings.T1_blood = 2.29; % T1 of arterial blood(s) ??
    settings.T2_blood = 68e-3; % T2 of blood(s) ??
    settings.T1_heart = 1.925; % T1 of heart(s) ??
    settings.T2_heart = 50e-3; % T2 of heart(s) ??
   
    mask = zeros(size(slice));
    % Define mask regions
    if settings.whole_body_mask == 0
        mask_heart = zeros(size(slice)); % pre-allocate heart mask
        mask_blood = zeros(size(slice)); % pre-allocate T2 value array
        for indi = 1:size(slice,1)
            for indj = 1:size(slice,2)
                if slice(indi,indj) > 0.85 && slice(indi,indj) < 0.95 % Assume Myocardium
                    mask_heart(indi,indj) = 1;
                elseif slice(indi,indj) > 1.3 && slice(indi,indj) < 1.4
                    mask_blood(indi,indj) = 1;
                end
            end
        end
        clear i j
        slice_T1 = mask_heart.*T1_heart + mask_blood.*T1_blood;
        slice_T2 = mask_heart.*T2_heart + mask_blood.*T2_blood;
        mask = mask_heart + mask_blood;
    elseif settings.Whole_body_mask == 1
        slice_T1 = 2.*slice; % Define slice T1 values as 2 * conductivity (very rough)
        slice_T1(round(slice,5,'significant') == round(slice(86,49),5,'significant')) =0; % Remove Lungs
        slice_T1 = slice_T1.*SyntheticDuke.mask(:,:,Syn_Slice); % Remove Tx/Rx
        slice_T2 = zeros(size(slice_T1));
        slice_T2(slice_T1 ~= 0) = 50e-3; % Give all values same T2
        mask(slice_T1 ~= 0) = 1;
    end
    if mTx == 0
        settings.Modes = 'NA';
        % Set up B1 Tx Field
        B1Tx = SyntheticDuke.B1Tx; % Transmit field 139 x 178 x 124 slices x 8 channels
        B1Tx_sumv = sum(bsxfun(@times,B1Tx(:,:,Syn_Slice,:),exp(1i*angle(conj(B1Tx(76,100,Syn_Slice,:))))),4); % Aligned fields to give max B1 field in centre of heart voxel (i=73,j=97) 53 103
        %B1Tx_sumv =sum(bsxfun(@times,B1Tx(:,:,Syn_Slice,:),permute(exp(1i.*[0,45,90,135,180,225,270,315]),[1,3,4,2])),4); % Sum of channels for Slice 60
        B1Tx_v = B1Tx_sumv/sqrt(50); % B1Tx in (T/V)
        B1Tx = B1Tx_v.*Ref_Voltage.*1e4; % Complex Magnitude of B1Tx in (G)
        Tx_FA_map = abs(B1Tx).*1e-4.*Gamma.*1e-3.*360; % Convert B1 map to FA

        imagesc(abs(Tx_FA_map(1:114,:))) % Plot B1Tx FA maps for each mode
        title('Simulated Flip Angle Map for a Synthetic Body Model')
        axis image off
        cb = colorbar;
        cb.Label.String = ['Flip angle, (',char(176),')'];
    elseif mTx == 1
        % Calculate Encoding Matrix
        [Enc_Mat,W_Mat] = Calc_Enc_Mat(settings.Enc_Scheme,settings.Modes);
        
        % Set up B1 Tx Field
        B1Tx = SyntheticDuke.B1Tx; % 139 x 178 x 124 slices x 8 channels
        B1Tx_modes = zeros(size(B1Tx,1),size(B1Tx,2),settings.Modes);
        Tx_FA_map = zeros(size(B1Tx_modes));
        for mode = 1:settings.Modes
            B1Tx_modes(:,:,mode) = sum(bsxfun(@times,B1Tx(:,:,Syn_Slice,:),permute(Enc_Mat(mode,:),[1,3,4,2])),4).*Ref_Voltage/sqrt(50); % Sum of channels for Slice 60 in (T)
            Tx_FA_map(:,:,mode) = (B1Tx_modes(:,:,mode)).*Gamma.*1e-3.*360; % Convert B1 map to FA (degrees) (this is complex)
        end
        imagesc(imtile(abs(Tx_FA_map),'Gridsize',[1 settings.Modes])) % Plot B1Tx FA maps for each mode
        title('Simulated Transmit Modes for Synthetic Body Model')
        axis image
        set(gca,'YTick',[]); set(gca,'XTick',[]);
        xlabel('Transmit Mode')
        xticks(((1:settings.Modes)).*size(Tx_FA_map,2) - size(Tx_FA_map,2)/2);
        xticklabels(1:settings.Modes);
        
        Tx_Channel_FA_map = squeeze(abs(B1Tx(:,:,Syn_Slice,:)).*Ref_Voltage.*Gamma.*1e-3.*360)/sqrt(50); % B1Tx in (T/V); % Convert B1 map to FA (degrees)
        %save('Tx_Channel_FA_map.mat','Tx_Channel_FA_map');
    end
    B1Rx = squeeze(SyntheticDuke.B1Rx(:,:,Syn_Slice,:)); % Receive field 139 x 178 x 124 slices x 8 channels
    
    % Don't re-simulate if input parameters unchanged
    %if
        % Function below handles simulating EPG image train and does analysis
        [results] = run_sequence_simulations(settings);
    %else
    %    disp('Input parameters unchanged - results not re-simulated.')
    %end
    Tx_FA_map = abs(Tx_FA_map); % now simulation is complete, absolute tx map applied (was complex)
    if  mTx == 1
        [outputArg1,outputArg2] = Pixelwise_Unencoding(Enc_Mat,W_Mat,settings.Modes)
    end


end
end



