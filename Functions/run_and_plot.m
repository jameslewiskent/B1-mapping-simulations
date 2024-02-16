function [results,settings] = run_and_plot(settings,plot_settings,results)
% This is the core function which controls plotting of the data. Hands off
% to other functions to simulate data if necessary- e.g. if not already simulated with identical parameters.

% ---------------------------------------------------------------------- %
% ----------- Following code handles simulation and plotting ----------- %
% ---------------------------------------------------------------------- %
settings.filename = ['Results_',char(datetime('now','TimeZone','local','Format','d-MMM-y-HH-mm')),'.mat'];
settings.savefilename = settings.filename; % In case of additional loop
settings.Scheme = lower(settings.Scheme); % change all scheme names to lowercase to prevent capitalisation from re-running otherwise identical simulations

if isempty(find(~isnan(settings.Noise), 1))
plot_settings.Dynamic_Range_Axis = 0;
warning('Dynamic range plot requested for data where no noise was simulated. The dynamic range calculated is not valid without noise present. Turning off dynamic range plotting.')
end

if settings.UseSyntheticData == 0
    Schemes = {'sattfl','sandwich','sa2rage','afi','dream'};
    if strcmpi(settings.Scheme,'ALL')
        settings = repmat(settings,1,length(Schemes));
        results = repmat(results,1,length(Schemes));
        for Scheme_n = 1:length(Schemes)
            disp([upper(Schemes{Scheme_n}),' (',num2str(Scheme_n),'/',num2str(length(Schemes)),')'])
            settings(Scheme_n).Scheme = Schemes{Scheme_n};
            settings(Scheme_n).MSor3D = 'default';
            [results(1,Scheme_n),settings(1,Scheme_n)] = run_sequence_simulations(settings(1,Scheme_n),results(1,Scheme_n),plot_settings); % Simulate and process data
        end
        
        places = [1,3,5,14,16];
        figure('color','w','units','normalized','position',[0.1,0.08,0.58,0.74],'Name','T1 Plot'); tiledlayout(4,6,'padding','none','tilespacing','compact');
        for Scheme_n = 1:length(Schemes)
            nexttile(places(Scheme_n),[2,2]); plot_T1_fig(results(1,Scheme_n),settings(1,Scheme_n),plot_settings); % T1 Plot
        end
        
        figure('color','w','units','normalized','position',[0.1,0.08,0.58,0.74],'Name','B0 Plot'); tiledlayout(4,6,'padding','none','tilespacing','compact');
        for Scheme_n = 1:length(Schemes)
            nexttile(places(Scheme_n),[2,2]); plot_B0_fig(results(1,Scheme_n),settings(1,Scheme_n),plot_settings); % B0 Plot
        end
        
        if size(settings(1,1).Velocities,2) > 1
            figure('color','w','units','normalized','position',[0.1,0.08,0.58,0.74],'Name','Flow Plot'); tiledlayout(4,6,'padding','none','tilespacing','compact');
            for Scheme_n = 1:length(Schemes)
                nexttile(places(Scheme_n),[2,2]); plot_Flow_fig(results(1,Scheme_n),settings(1,Scheme_n),plot_settings);
            end
        end
        
        if size(settings(1,1).Diff_coeffs,2) > 1
            figure('color','w','units','normalized','position',[0.1,0.08,0.58,0.74],'Name','Diff Plot'); tiledlayout(4,6,'padding','none','tilespacing','compact');
            for Scheme_n = 1:length(Schemes)
                nexttile(places(Scheme_n),[2,2]); plot_Diff_fig(results(1,Scheme_n),settings(1,Scheme_n),plot_settings);
            end
        end
        
        plot_SAR_fig(results,settings); % SAR Plot
        if settings(1,1).Calc_FWHM ~= 0
        plot_FWHM(results,settings,plot_settings);
        end
    else
        
        if isfield(settings,'LoopValues') % Can add other loops here if you wish
            for Additional_Loop_Counter = 1:length(settings.LoopValues)
                settings.Additional_Loop_Counter = Additional_Loop_Counter;
                [results,settings] = run_sequence_simulations(settings,results,plot_settings); % Simulate and process data
            end
        else
            [results,settings] = run_sequence_simulations(settings,results,plot_settings); % Simulate and process data
        end
        if isfield(results,'Dynamic_Range')
        disp(['Dynamic Range: ',num2str(results.Dynamic_Range,'%4.1f')]);
        end
        
        if settings.HR_TR == 1
            figure('color','w');
            nexttile();
            histogram(results.seq_TRs.*1e3,100)
            xlim([0 max(results.seq_TRs.*1e3,[],'all')])
            xticks([0:250:max(results.seq_TRs.*1e3,[],'all')])
            xlabel('TR [ms]'); ylabel('Frequency');
            title('HRV Histogram')
            hold on;
            histogram(abs(diff(results.seq_TRs,1)).*1e3,100)
            xlabel('Time [ms]'); ylabel('Frequency');
            axis square
            legend({'TR','HRV'})
            
            
        end
        
        figure('color','w','Name','T1 Plot'); plot_T1_fig(results,settings,plot_settings); % T1 Plot
        figure('color','w','Name','B0 Plot'); plot_B0_fig(results,settings,plot_settings); % B0 Plot
        if size(settings.Velocities,2) > 1
            figure('color','w','Name','Flow Plot'); plot_Flow_fig(results,settings,plot_settings); % Flow Plot
        end
        if size(settings.Diff_coeffs,2) > 1
            figure('color','w','Name','Diff Plot'); plot_Diff_fig(results,settings,plot_settings); % Flow Plot
        end
        figure('color','w','Name','Magnetisation Plot'); plot_mz(results,settings); % magnetisation plot
        figure('color','w','Name','RF Pulse Plot'); plot_rf_pulses(settings)
        
        if settings.Calc_FWHM ~= 0
        plot_FWHM(results,settings,plot_settings);
        end
    end
else
    % ---------------------------------------------------------------------- %
    % -----------             Simulating Synthetic Data          ----------- %
    % ---------------------------------------------------------------------- %
    
    % Read in and generate slice data if not already existing
    SyntheticDuke = load('Data\SyntheticDuke.mat'); % Reads in if current folder is Masterscript
    % Assign T1, T2 values to synthetic data, mask out non-heart regions (temp - to speed up dev)
    
    slice = SyntheticDuke.sigma(:,:,settings.Syn_Slice); % For now only consider one slice of conductivity data to mask
    
    % Define relaxation constants for blood & heart
    settings.T1_blood = 2.29; % T1 of arterial blood(s) ??
    settings.T2_blood = 68e-3; % T2 of blood(s) ??
    settings.T1_heart = 1.925; % T1 of heart(s) ??
    settings.T2_heart = 50e-3; % T2 of heart(s) ??
    
    mask = zeros(size(slice));
    % Define mask regions
    if settings.Whole_body_mask == 0
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
        slice_T1 = slice_T1.*SyntheticDuke.mask(:,:,settings.Syn_Slice); % Remove Tx/Rx
        slice_T2 = zeros(size(slice_T1));
        slice_T2(slice_T1 ~= 0) = 50e-3; % Give all values same T2
        mask(slice_T1 ~= 0) = 1;
    end

    B1Rx = squeeze(SyntheticDuke.B1Rx(:,:,settings.Syn_Slice,:)); % Receive field 139 x 178 x 124 slices x 8 channels
    
    % Don't re-simulate if input parameters unchanged
    %if
    % Function below handles simulating EPG image train and does analysis
    [results] = run_sequence_simulations(settings,results,plot_settings);
    %else
    %    disp('Input parameters unchanged - results not re-simulated.')
    %end
    Tx_FA_map = abs(Tx_FA_map); % now simulation is complete, absolute tx map applied (was complex)
    if  mTx == 1
        [outputArg1,outputArg2] = Pixelwise_Unencoding(Enc_Mat,W_Mat,settings.Modes)
    end
    
end
end

