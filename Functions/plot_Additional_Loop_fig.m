function plot_Additional_Loop_fig(results,settings,plot_settings)
% Additional Loop Plots
Noise_n = find(~isnan(settings.Noise)); % First non-NaN (non-zero noise)
if isempty(Noise_n)
    Noise_n = 1;
else
    Noise_n = Noise_n(1);
end

% Plot lookup table for each TR
if strcmp(settings.LoopFieldName,'TR')
    % Load in various files from the loop
    T1_n = 1:length(settings.T1s);
    B0_n = 1;
    Flow_n = 1;
    Diff_n = 1;
    for Additional_Loop_Counter = 1:length(settings.LoopValues)
        settings.(settings.LoopFieldName) = settings.LoopValues(Additional_Loop_Counter,:);
        [already_ran,filename] = check_if_already_ran(settings);
        if ~already_ran
            error('Cannot find results to plot!')
            clear filename
        end
        
        load([settings.filepath,filesep,filename],'results');
        FAs(Additional_Loop_Counter,:,:,:) = squeeze(mean(results.Measured_FA(1,1,1,:,B0_n,T1_n,Flow_n,Diff_n,Noise_n,:),[10]));
        DR(Additional_Loop_Counter) = results.Dynamic_Range;
        DR_Values(Additional_Loop_Counter,:) = results.DR_Values;
    end
    DR_ratio = DR./settings.LoopValues';
    figure('color','w','Units','normalized','position',[0.2,0.3,0.575,0.486]);
    plot(settings.LoopValues,DR,'rx'); hold on
    xlim([0.5 2])
    ylim([1 10])
    %ylim([floor(min(DR_ratio)*10)/10 ceil(max(DR_ratio)*10)/10])
    xticks(settings.LoopValues(1):0.2:settings.LoopValues(end))
    ylabel('Dynamic Range Efficiency, [s^{-1}]'); xlabel('TR, [s]')
    grid on
    axis square
end


if strcmp(settings.LoopFieldName,'Matrix_Size')
    % Plot lookup table for each Matrix Size
    % Load in various files from the loop
    for Additional_Loop_Counter = 1:length(settings.LoopValues)
        lookup_filename = [settings.Scheme,'_lookup_table_Loop',num2str(Additional_Loop_Counter),'.mat'];
        load([settings.filepath,filesep,lookup_filename],'x_query','fx_interp');
        fx_interps(Additional_Loop_Counter,:) = fx_interp;
        x_queries(Additional_Loop_Counter,:) = x_query;
    end
end

if strcmp(settings.LoopFieldName,'Ejection_Fraction')
    Nominal_FA = settings.Dynamic_Range.*settings.nomPP_FA*(180/pi);
    lower_bound = Nominal_FA(35);
    upper_bound = Nominal_FA(160);
    T1_n = 1:length(settings.T1s);
    B0_n = 1;
    Flow_n = 1;
    Diff_n = 1;
    
    for Additional_Loop_Counter = 1:length(settings.LoopValues)
        settings.(settings.LoopFieldName) = settings.LoopValues(Additional_Loop_Counter,:);
        [already_ran,filename] = check_if_already_ran(settings);
        if ~already_ran
            error('Cannot find results to plot!')
            clear filename
        end
        
        load([settings.filepath,filesep,filename],'results');
        FAs(Additional_Loop_Counter,:,:,:) = squeeze(mean(results.Measured_FA(1,1,1,:,B0_n,T1_n,Flow_n,Diff_n,Noise_n,:),10));
    end
    
    Diff_FAs = (FAs(:,:,:) - Nominal_FA(1,:))./Nominal_FA(1,:); % Difference to EF = 0
    
    Mean_Diff_FAs = 100.*squeeze(mean(Diff_FAs,3));
    figure('color','w'); tiledlayout(1,2,'padding','none','tilespacing','compact');
    nexttile;
    imagesc(Nominal_FA,settings.LoopValues.*100,Mean_Diff_FAs,[-10 10]);
    set(gca,'YDir','normal')
    ylabel('Ejection Fraction, [%]');
    xlabel(['Nominal FA, [',char(176),']']);
    xlim([0 200]);
    axis square
    set(gca,'Ydir','normal')
    colormap(bluewhitered)
    cb = colorbar; cb.Label.String = 'Percent Error, [%]';
    xline(lower_bound,'linewidth',2,'linestyle','--');
    xline(upper_bound,'linewidth',2,'linestyle','--');
    title('a) Mean T1');
    
    
    Mean_Diff_FAs = 100.*squeeze(mean(Diff_FAs(:,lower_bound:upper_bound,:),2));
    nexttile();
    imagesc(settings.T1s,settings.LoopValues.*100,Mean_Diff_FAs,[-10 10]);
    set(gca,'YDir','normal')
    ylabel('Ejection Fraction, [%]');
    xlabel('T_1, [s]');
    axis square
    set(gca,'Ydir','normal')
    colormap(bluewhitered)
    cb = colorbar; cb.Label.String = 'Percent Error, [%]';
    title('b) Mean DR');
end

if strcmp(settings.LoopFieldName,'Coil_Cycle')
    
    for Additional_Loop_Counter = 1:size(settings.LoopValues,1)
        for Additional_Loop_Counter2 = 1:size(settings.LoopValues2,1)
            settings.(settings.LoopFieldName) = settings.LoopValues(Additional_Loop_Counter,:);
            settings.(settings.LoopFieldName2) = settings.LoopValues2(Additional_Loop_Counter2,:);
            
            if strcmpi(settings.Scheme,'GRE') && ~strcmpi(settings.LoopFieldName2,'Coil_Cycle_Order')
                if settings.Coil_Cycle == 0
                    settings.Coil_Cycle_Order = 1:8;
                elseif settings.Coil_Cycle == 1
                    settings.Coil_Cycle_Order = [1,4,7,2,5,8,3,6];
                end
            end
            
            % This is here incase FA is changed, as RF pulse will need updating too
            settings.RF_Pulse = Get_RF_Pulse(settings.nom_FA,settings.RF_Type,settings.RF_Time,settings.RF_TBP,settings.Ref_Voltage);
            
            % If Matrix size changed, these needs updating
            settings.Scan_Size(1) = round(settings.PE1_Resolution*settings.PE1_Partial_Fourier*settings.Matrix_Size(1));
            settings.Scan_Size(2) = round(settings.PE2_Resolution*settings.PE2_Partial_Fourier*settings.Matrix_Size(2));
            settings = Calc_Slice_Positions(settings); % Calculate Slice positions
            settings = Calc_Segment_Sizes(settings);
            settings.prep_spoils = 1.* ones(1,settings.N_TRs); % Number of unit gradients to move through after pre-pulse (spoiling- need enough for optional dummy scans too)
            settings.train_spoils = 1.* ones(1,ceil(settings.Scan_Size(1)/settings.Segment_Factor));
            
            if settings.UseSyntheticData == 1
                [settings.Dynamic_Range,settings.Tx_FA_map,settings.Enc_Mat] = Calc_Tx(max(settings.RF_Pulse),settings); % Now pulse voltages are known, we can calculate the associated transmit field
            end
            
            [already_ran,filename] = check_if_already_ran(settings);
            if ~already_ran
                error('Cannot find results to plot!')
                clear filename
            end
            clear results;
            load([settings.filepath,filesep,filename],'results');
            
            Rel_Image_Maps(:,:,:,:,:,:,:,:,:,Additional_Loop_Counter,Additional_Loop_Counter2) = mean(results.Image_Maps(:,:,:,1,:,:,:,:,:,:),10);
            
            % Calculate ground truth relative maps
            GT_Rel_Image_Maps(:,:,:,Additional_Loop_Counter,Additional_Loop_Counter2) = settings.Tx_FA_map./sum(abs(settings.Tx_FA_map),3);
        end
    end
    
    % Check whether maps both have values in same voxels
    x = GT_Rel_Image_Maps;
    x(x ~= 0) = 1;
    x = ~any(~x,3:ndims(x)); % Must check extra dimensions for simulated maps, then collapse
    y = Rel_Image_Maps;
    y(y ~= 0) = 1;
    y = ~any(~y,3:ndims(y)); % Must check extra dimensions for simulated maps, then collapse
    diff_mask = logical(x-y);
    % Set to NaN any voxels which aren't in both the ground truth and simulated maps
    GT_Rel_Image_Maps(repmat(diff_mask,[1 1 size(GT_Rel_Image_Maps,3:ndims(GT_Rel_Image_Maps))])) = NaN;
    Rel_Image_Maps(repmat(diff_mask,[1 1 size(Rel_Image_Maps,3:ndims(Rel_Image_Maps))])) = NaN;
    
    % Create long arrays & remove NaNs
    for Tx_n = 1:8
        for Additional_Loop_Counter = 1:size(settings.LoopValues,1)
            for Additional_Loop_Counter2 = 1:size(settings.LoopValues2,1)
                for B0_n = 1:length(settings.B0_Range_Hz)
                    for T1_n = 1:length(settings.T1s)
                        for Flow_n = 1:length(settings.Velocities)
                            for Diff_n = 1:length(settings.Diff_coeffs)
                                for Noise_n = 1:length(settings.Noise)
                                    % Remove NaNs from ground truth
                                    Masked_GT_Rel_Image_Maps = settings.Synthetic_Mask.*abs(GT_Rel_Image_Maps(:,:,Tx_n,Additional_Loop_Counter,Additional_Loop_Counter2));
                                    Masked_GT_Rel_Image_Maps(Masked_GT_Rel_Image_Maps == 0) = NaN;
                                    GT_Rel_Long_Maps(:,Tx_n,Additional_Loop_Counter,Additional_Loop_Counter2) = Masked_GT_Rel_Image_Maps(~isnan(Masked_GT_Rel_Image_Maps));
                                    % Remove NaNs from results
                                    Masked_Rel_Image_Maps = settings.Synthetic_Mask.*abs(Rel_Image_Maps(:,:,Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2));
                                    Masked_Rel_Image_Maps(Masked_Rel_Image_Maps == 0) = NaN;
                                    Rel_Long_Maps(:,Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2) = Masked_Rel_Image_Maps(~isnan(Masked_Rel_Image_Maps));
                                    % Calculate NRMSE
                                    NRMSE_Values(Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2) = rmse(Rel_Long_Maps(:,Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2),GT_Rel_Long_Maps(:,Tx_n,Additional_Loop_Counter,Additional_Loop_Counter2) ,'norm');
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    % Values to plot
    T1_n = 1;
    B0_n = 1;
    Flow_n = 1;
    Diff_n = 1;
    
    % Maps
    figure('color','w','Units','Normalized','Position',[0.197395833333333,0.142592592592593,0.571354166666667,0.713888888888888])
    tiledlayout('flow','tilespacing','none','padding','none');
    for Additional_Loop_Counter = 1:size(settings.LoopValues,1)
        for Additional_Loop_Counter2 = 1:size(settings.LoopValues2,1)
            nexttile;    imagesc(imtile(abs(Rel_Image_Maps(:,:,:,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2)),'GridSize', [2 4]),[0 1])
            axis image off
            set(gca,'YDir','normal');
        end
    end
    set(gca,'YDir','normal');
    cb = colorbar;
    cb.Layout.Tile = 'East'; % <-- place legend east of tiles
    cb.Label.String = 'Relative B_1^+ [a.u.]';
    cb.FontSize = 16;
    
    % Difference
    figure('color','w','Units','Normalized','Position',[0.197395833333333,0.142592592592593,0.571354166666667,0.713888888888888])
    tiledlayout('flow','tilespacing','none','padding','none');
    for Additional_Loop_Counter = 1:size(settings.LoopValues,1)
        for Additional_Loop_Counter2 = 1:size(settings.LoopValues2,1)
            data_to_plot = 100.*settings.Body_Mask.*(abs(Rel_Image_Maps(:,:,:,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2)) - abs(GT_Rel_Image_Maps(:,:,:,Additional_Loop_Counter,Additional_Loop_Counter2)))./abs(GT_Rel_Image_Maps(:,:,:,Additional_Loop_Counter,Additional_Loop_Counter2));
            data_to_plot(isnan(data_to_plot)) = 0;
            nexttile;    imagesc(imtile(data_to_plot,'GridSize', [2 4]),[-100 100])
            axis image off
            set(gca,'YDir','normal');
        end
    end
    colormap(bluewhitered)
    cb = colorbar;
    cb.Layout.Tile = 'East'; % <-- place legend east of tiles
    cb.Label.String = 'Relative Difference [%]';
    cb.FontSize = 16;
    
    if ~isfield(settings,'LoopFieldName2') || strcmp(settings.LoopFieldName2,'Coil_Cycle_Order')
        % Correlation plots
        xlimmin = 0; xlimmax = 1;
        ylimmin = 0; ylimmax = 1;
        Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
        fontsize = 10; Font = 'calibri'; Axis_fontsize = 11;
        
        figure('color','w','Units','normalized','Position',[0.317057291666667,0.047222222222222,0.461067708333333,0.857407407407407])
        cmap = jet(8);
        tiledlayout('flow','tilespacing','compact','padding','none');
        n = 1;
        for Additional_Loop_Counter = 1:size(settings.LoopValues,1)
            for Additional_Loop_Counter2 = 1:size(settings.LoopValues2,1)
                
                if strcmpi(settings.LoopFieldName,'Coil_Cycle') && settings.LoopValues(Additional_Loop_Counter) == 0
                    titlepart1 = 'Channelwise ';
                elseif strcmpi(settings.LoopFieldName,'Coil_Cycle') && settings.LoopValues(Additional_Loop_Counter) == 1
                    titlepart1 = 'Coil-cycled ';
                else
                    titlepart1 = ' ';
                end
                if strcmpi(settings.LoopFieldName2,'Coil_Cycle_Order') && all(settings.LoopValues2(Additional_Loop_Counter2,:) == 1:8)
                    titlepart2 = ' ';
                elseif strcmpi(settings.LoopFieldName2,'Coil_Cycle_Order') && all(settings.LoopValues2(Additional_Loop_Counter2,:) == [1,4,7,2,5,8,3,6])
                    titlepart2 = 'Interleaved ';
                else
                    titlepart2 = ' ';
                end
                fulltitle = [titlepart1,titlepart2];
                
                nexttile; textx = xlimmin + 0.05; texty = ylimmax - 0.05;
                for Tx_n = 1:8
                    plot(abs(Rel_Long_Maps(:,Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2)),abs(GT_Rel_Long_Maps(:,Tx_n,Additional_Loop_Counter,Additional_Loop_Counter2)),'color',cmap(Tx_n,:),'Marker','.','Linestyle','none','Markersize',2,'HandleVisibility','off'); hold on
                end
                plot([xlimmin,xlimmax],[ylimmin,ylimmax],'k','HandleVisibility','off')
                
                Long_Maps = reshape(abs(Rel_Long_Maps(:,:,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2)),[],size(settings.LoopValues,1),size(settings.LoopValues2,1));
                Long_GT_Maps = reshape(abs(GT_Rel_Long_Maps(:,:,Additional_Loop_Counter,Additional_Loop_Counter2)),[],size(settings.LoopValues,1),size(settings.LoopValues2,1));
                
                P = polyfit(Long_Maps(:,Additional_Loop_Counter,Additional_Loop_Counter2),Long_GT_Maps(:,Additional_Loop_Counter,Additional_Loop_Counter2),1);
                text(textx,texty,['y = ',num2strpn(P(1),'%4.2f'),'x ',num2strpn(P(2),'%4.2f')],'Fontsize',fontsize); texty = texty - 0.05;
                [~,~,~,~,stats] = regress(Long_Maps(:,Additional_Loop_Counter,Additional_Loop_Counter2),Long_GT_Maps(:,Additional_Loop_Counter,Additional_Loop_Counter2));
                text(textx,texty,['r^2 = ',num2str(stats(1),'%4.2f')],'Fontsize',fontsize);  texty = texty - 0.05;
                text(textx,texty,['n = ',num2str(size(Long_Maps,1))],'Fontsize',fontsize);  texty = texty - 0.05;
                text(textx,texty,['RMSE = ',num2str(mean(NRMSE_Values(:,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2),1),'%4.2f')],'Fontsize',fontsize);
                
                title([Alphabet(n),') ',fulltitle],'FontName',Font,'Fontsize',16)
                xlabel('Ground-Truth [a.u.]','FontName',Font,'FontSize',Axis_fontsize,'FontWeight','normal')
                ylabel('Relative B_1^+ [a.u.]','FontName',Font,'FontSize',Axis_fontsize,'FontWeight','normal')
                box on
                axis square
                xlim([xlimmin xlimmax])
                ylim([ylimmin ylimmax])
                
                if Additional_Loop_Counter == 1 && Additional_Loop_Counter2 == 1
                    for Tx_n = 1:8
                        plot(-10,-10,'color',cmap(Tx_n,:),'Marker','.','LineStyle','none','MarkerSize',20)
                    end
                    leg = legend('1','2','3','4','5','6','7','8','Location','southeast');
                    title(leg,'Tx #');
                end
                n = n + 1;
            end
        end
    end
    
    if strcmp(settings.LoopFieldName2,'nom_FA')
        figure('color','w','Units','Normalized','Position',[0.317057291666667,0.188888888888889,0.416276041666667,0.691898148148148])
        cmap = jet(8); Markers = {'x','o'};
        for Additional_Loop_Counter = 1:size(settings.LoopValues,1)
            for Additional_Loop_Counter2 = 1:size(settings.LoopValues2,1)
                for Tx_n = 1:8
                    plot(settings.LoopValues2(Additional_Loop_Counter2).*180/pi,100.*NRMSE_Values(Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2),'color',cmap(Tx_n,:),'Marker',Markers{Additional_Loop_Counter},'MarkerSize',10,'Linewidth',2,'HandleVisibility','off'); hold on
                end
            end
        end
        plot(-10,-10,'color','k','Marker',Markers{1},'LineStyle','none','MarkerSize',10,'Linewidth',2)
        plot(-10,-10,'color','k','Marker',Markers{2},'LineStyle','none','MarkerSize',10,'Linewidth',2)
        legend('Channel-wise','Coil-cycled','Location','NorthWest');
        xlim([0 15]); ylim([0 100]);
        grid on
        axis square
        ylabel('NRMSE, [%]')
        xlabel(['Nominal Flip Angle, [',char(176),']']);
    end
    
    if size(settings.T1s,2) > 1
        figure('color','w','Units','Normalized','Position',[0.317057291666667,0.188888888888889,0.416276041666667,0.691898148148148])
        cmap = jet(8); Markers = {'x','o'};
        for Additional_Loop_Counter = 1:size(settings.LoopValues,1)
            for T1_n = 1:size(settings.T1s,2)
                for Tx_n = 1:8
                    plot(settings.T1s(T1_n),100.*NRMSE_Values(Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2),'color',cmap(Tx_n,:),'Marker',Markers{Additional_Loop_Counter},'MarkerSize',10,'Linewidth',2,'HandleVisibility','off'); hold on
                end
            end
        end
        plot(-10,-10,'color','k','Marker',Markers{1},'LineStyle','none','MarkerSize',10,'Linewidth',2)
        plot(-10,-10,'color','k','Marker',Markers{2},'LineStyle','none','MarkerSize',10,'Linewidth',2)
        legend('Channel-wise','Coil-cycled','Location','NorthWest');
        xlim([0 4]); ylim([0 100]);
        grid on
        axis square
        ylabel('NRMSE, [%]')
        xlabel('T_1, [s]');
    end
    
    if strcmp(settings.LoopFieldName2,'Matrix_Size')
        figure('color','w','Units','Normalized','Position',[0.317057291666667,0.188888888888889,0.416276041666667,0.691898148148148])
        cmap = jet(8); Markers = {'x','o'};
        for Additional_Loop_Counter = 1:size(settings.LoopValues,1)
            for Additional_Loop_Counter2 = 1:size(settings.LoopValues2,1)
                for Tx_n = 1:8
                    plot(settings.LoopValues2(Additional_Loop_Counter2),100.*NRMSE_Values(Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2),'color',cmap(Tx_n,:),'Marker',Markers{Additional_Loop_Counter},'MarkerSize',10,'Linewidth',2,'HandleVisibility','off'); hold on
                end
            end
        end
        plot(-10,-10,'color','k','Marker',Markers{1},'LineStyle','none','MarkerSize',10,'Linewidth',2)
        plot(-10,-10,'color','k','Marker',Markers{2},'LineStyle','none','MarkerSize',10,'Linewidth',2)
        legend('Channel-wise','Coil-cycled','Location','NorthWest');
        xlim([0 72]); ylim([0 100]);
        grid on
        axis square
        ylabel('NRMSE, [%]')
        xlabel('N_{PE}, [s]');
    end
    
    
end


end