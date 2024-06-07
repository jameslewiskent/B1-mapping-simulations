function [results,settings] = run_and_plot(settings,plot_settings,results)
% This is the core function which controls plotting of the data. Hands off
% to other functions to simulate data if necessary- e.g. if not already simulated with identical parameters.

% ---------------------------------------------------------------------- %
% ----------- Following code handles simulation and plotting ----------- %
% ---------------------------------------------------------------------- %
settings.Scheme = lower(settings.Scheme); % change all scheme names to lowercase to prevent capitalisation from re-running otherwise identical simulations

if isempty(find(~isnan(settings.Noise), 1))
    plot_settings.Dynamic_Range_Axis = 0;
    warning('Dynamic range plot requested for data where no noise was simulated. The dynamic range calculated is not valid without noise present. Turning off dynamic range plotting.')
end

if ~isfield(settings,'LoopValues2')
    settings.LoopFieldName2 = 'NaN'; 
    settings.LoopValues2 = NaN;
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
            for Additional_Loop_Counter = 1:size(settings.LoopValues,1)
                settings.Additional_Loop_Counter = Additional_Loop_Counter;
                
                for Additional_Loop_Counter2 = 1:size(settings.LoopValues2,1)
                    settings.Additional_Loop_Counter2 = Additional_Loop_Counter2;
                    
                    disp(['Innerloop: ',num2str(Additional_Loop_Counter),', Outerloop: ',num2str(Additional_Loop_Counter2)])
                    [results,settings] = run_sequence_simulations(settings,results,plot_settings); % Simulate and process data
                end
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
            xticks(0:250:max(results.seq_TRs.*1e3,[],'all'))
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
    
    [settings.Body_Mask, settings.Heart_Mask, settings.Synthetic_T1s] = Create_Synthetic_Masks(settings.Syn_Slice,settings);
    
    % [settings] = Choose_Mag_Track_Locations(settings)
    
    % Choose to use whole body or heart
    if settings.Whole_body_mask == 0
        settings.Synthetic_Mask = settings.Heart_Mask;
    elseif settings.Whole_body_mask == 1
        settings.Synthetic_Mask = settings.Body_Mask;
    end
    
    if settings.verbose == 1
        % Plot synthetic T1 figure
        figure('color','w'); tiledlayout('flow','Tilespacing','none','Padding','none'); nexttile;
        imagesc(settings.Synthetic_T1s); hold on
        set(gca,'Ydir','normal')
        colormap(jet)
        axis image off
        cb = colorbar;
        cb.Label.String = 'T_1 [s]';
        
        if size(settings.Mag_Track_SynInd,1) > 1
            % Plot positions of magnetisation tracking onto T1 figure
            %cmap = hsv(size(settings.Mag_Track_SynInd,1));
            cmap = repmat([1 0 0],size(settings.Mag_Track_SynInd,1),1);
            for location_n = 1:size(settings.Mag_Track_SynInd,1) 
                plot(settings.Mag_Track_SynInd(location_n,2),settings.Mag_Track_SynInd(location_n,1),'color',cmap(location_n,:),'Marker','x','MarkerSize',10,'LineWidth',3); hold on
            end
            hold off
        end
    end
    
    settings.Long_Synthetic_T1s = reshape(settings.Synthetic_T1s,[],1);
    settings.Long_Synthetic_Mask = reshape(settings.Synthetic_Mask,[],1);
    %B1Rx = squeeze(SyntheticDuke.B1Rx(:,:,settings.Syn_Slice,:)); % Receive field 139 x 178 x 124 slices x 8 channels
    
    mz_h = figure('color','w','Name','Magnetisation Plot'); tiledlayout(mz_h,'flow','tilespacing','compact','padding','none');
    tx_h = figure('color','w','Name','Magnetisation Plot'); tiledlayout(tx_h,'flow','tilespacing','compact','padding','none');
    if isfield(settings,'LoopValues') % Can add other loops here if you wish
        for Additional_Loop_Counter = 1:size(settings.LoopValues,1)
            settings.Additional_Loop_Counter = Additional_Loop_Counter;
            
            for Additional_Loop_Counter2 = 1:size(settings.LoopValues2,1)
                settings.Additional_Loop_Counter2 = Additional_Loop_Counter2;
                
                disp(['Innerloop: ',num2str(Additional_Loop_Counter),', Outerloop: ',num2str(Additional_Loop_Counter2)])
                [results,settings] = run_sequence_simulations(settings,results,plot_settings); % Simulate and process data
                
                figure(mz_h); nexttile; plot_mz(results,settings); % magnetisation plot
                figure(tx_h); nexttile; plot_tx(results,settings); % transmit field plot
            end
        end
        plot_Additional_Loop_fig(results,settings,plot_settings)
    else
        % Function below handles simulating EPG image train and does analysis
        [results] = run_sequence_simulations(settings,results,plot_settings);
        
        % Plot image maps
        Noise_n = 2;
        figure('color','w');
        imagesc(imtile(abs(results.Image_Maps(:,:,:,1,1,1,1,1,Noise_n,1)),'Gridsize',[1 settings.Modes])) % Plot B1Tx FA maps for each mode
        set(gca,'Ydir','normal')
        axis image
        set(gca,'YTick',[]); set(gca,'XTick',[]);
        colormap(turbo)
        colorbar
        
        figure('color','w','Name','Magnetisation Plot'); plot_mz(results,settings); % magnetisation plot
    end
end
end

