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
        if strcmpi(settings.LoopFieldName,'Coil_Cycle') && strcmpi(settings.LoopValues(Additional_Loop_Counter),'SS')
            titlepart1{Additional_Loop_Counter} = 'Channel-wise ';
        elseif strcmpi(settings.LoopFieldName,'Coil_Cycle') && strcmpi(settings.LoopValues(Additional_Loop_Counter),'SW')
            titlepart1{Additional_Loop_Counter} = 'Channel-wise ';
        elseif strcmpi(settings.LoopFieldName,'Coil_Cycle') && strcmpi(settings.LoopValues(Additional_Loop_Counter),'CC')
            titlepart1{Additional_Loop_Counter} = 'Coil-cycled ';
        else
            titlepart1{Additional_Loop_Counter} = ' ';
        end
    end
    
    for Additional_Loop_Counter = 1:size(settings.LoopValues,1)
        for Additional_Loop_Counter2 = 1:size(settings.LoopValues2,1)
            settings.(settings.LoopFieldName) = settings.LoopValues(Additional_Loop_Counter,:);
            settings.(settings.LoopFieldName2) = settings.LoopValues2(Additional_Loop_Counter2,:);
            
            % In case of default schemes
            if strcmpi(settings.Scheme,'GRE') && ~strcmpi(settings.LoopFieldName2,'Coil_Cycle_Order')
                if strcmpi(settings.Coil_Cycle,'SS') || strcmpi(settings.Coil_Cycle,'SW')
                    settings.Coil_Cycle_Order = 1:8;
                elseif strcmpi(settings.Coil_Cycle,'CC')
                    settings.Coil_Cycle_Order = [1,4,7,2,5,8,3,6];
                end
            end
            
            % This is here incase FA is changed, as RF pulse will need
            % updating too
            settings.RF_Pulse = Get_RF_Pulse(settings.nom_FA,settings.RF_Type,settings.RF_Time,settings.RF_TBP,settings.Ref_Voltage);
            
            % If Matrix size changed, these needs updating
            settings.Scan_Size(1) = round(settings.PE1_Resolution*settings.PE1_Partial_Fourier*settings.Matrix_Size(1));
            settings.Scan_Size(2) = round(settings.PE2_Resolution*settings.PE2_Partial_Fourier*settings.Matrix_Size(2));
            settings = Calc_Slice_Positions(settings); % Calculate Slice positions
            settings = Calc_Segment_Sizes(settings);
            settings.prep_spoils = 1.* ones(1,settings.N_TRs); % Number of unit gradients to move through after pre-pulse (spoiling- need enough for optional dummy scans too)
            settings.train_spoils = 1.* ones(1,ceil(settings.Scan_Size(1)/settings.Segment_Factor));
            settings = Check_Min_TR(settings); % Check minimum TR time (Not currently used by GRE sequence)
            
            if (strcmpi(settings.UseSyntheticData,'Duke') || strcmpi(settings.UseSyntheticData,'Phantom'))
                % These are the ground truth mode transmit maps
                [settings.Dynamic_Range,settings.Tx_FA_map,settings.Enc_Mat] = Calc_Tx(max(settings.RF_Pulse),settings); % Now pulse voltages are known, we can calculate the associated transmit field
                
                % These are the ground truth single-channel maps
                if settings.Modes > 1 && ~strcmpi(settings.Enc_Scheme,'Indiv')
                    [~,settings.GT_Tx_FA_map,~] = Calc_Tx(max(settings.RF_Pulse),settings,Calc_Enc_Mat('Indiv',settings.Modes)); % Now pulse voltages are known, we can calculate the associated transmit field
                else
                    settings.GT_Tx_FA_map = settings.Tx_FA_map;
                end
            end
            
            [already_ran,filename] = check_if_already_ran(settings);
            if ~already_ran
                error('Cannot find results to plot!')
                clear filename
            end
            clear results;
            load([settings.filepath,filesep,filename],'results');
            
            Rel_Images(:,:,:,:,:,:,:,:,:,Additional_Loop_Counter,Additional_Loop_Counter2) = reshape(permute(results.Max_Val_IT1,[4,2,3,1,5,6,7,8,9]),[size(settings.Tx_FA_map),1,size(results.Measured_FA,5:ndims(results.Measured_FA))]);
            
            Rel_Image_Maps(:,:,:,:,:,:,:,:,:,Additional_Loop_Counter,Additional_Loop_Counter2) = mean(results.Image_Maps(:,:,:,1,:,:,:,:,:,:),10);
            
            % Calculate ground truth relative maps
            GT_Rel_Image_Maps(:,:,:,Additional_Loop_Counter,Additional_Loop_Counter2) = settings.GT_Tx_FA_map./sum(abs(settings.GT_Tx_FA_map),3);
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
                                    
                                    if strcmpi(settings.UseSyntheticData,'Duke')
                                        % Remove NaNs from heart masked ground truth
                                        Heart_Masked_GT_Rel_Image_Maps = settings.Heart_Mask.*abs(GT_Rel_Image_Maps(:,:,Tx_n,Additional_Loop_Counter,Additional_Loop_Counter2));
                                        Heart_Masked_GT_Rel_Image_Maps(Heart_Masked_GT_Rel_Image_Maps == 0) = NaN;
                                        Heart_GT_Rel_Long_Maps(:,Tx_n,Additional_Loop_Counter,Additional_Loop_Counter2) = Heart_Masked_GT_Rel_Image_Maps(~isnan(Heart_Masked_GT_Rel_Image_Maps));
                                    end
                                    
                                    % Remove NaNs from results
                                    Masked_Rel_Image_Maps = settings.Synthetic_Mask.*abs(Rel_Image_Maps(:,:,Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2));
                                    Masked_Rel_Image_Maps(Masked_Rel_Image_Maps == 0) = NaN;
                                    Rel_Long_Maps(:,Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2) = Masked_Rel_Image_Maps(~isnan(Masked_Rel_Image_Maps));
                                    
                                    if strcmpi(settings.UseSyntheticData,'Duke')
                                        % Remove NaNs from heart masked results
                                        Heart_Masked_Rel_Image_Maps = settings.Heart_Mask.*abs(Rel_Image_Maps(:,:,Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2));
                                        Heart_Masked_Rel_Image_Maps(Heart_Masked_Rel_Image_Maps == 0) = NaN;
                                        Heart_Rel_Long_Maps(:,Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2) = Heart_Masked_Rel_Image_Maps(~isnan(Heart_Masked_Rel_Image_Maps));
                                    end
                                    
                                    % Calculate NRMSE
                                    NRMSE_Values(Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2) = rmse(Rel_Long_Maps(:,Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2),GT_Rel_Long_Maps(:,Tx_n,Additional_Loop_Counter,Additional_Loop_Counter2) ,'norm');
                                    
                                    if strcmpi(settings.UseSyntheticData,'Duke')
                                        % Calculate heart masked NRMSE
                                        Heart_NRMSE_Values(Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2) = rmse(Heart_Rel_Long_Maps(:,Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2),Heart_GT_Rel_Long_Maps(:,Tx_n,Additional_Loop_Counter,Additional_Loop_Counter2) ,'norm');
                                    end
                                    
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
    
    if ~isfield(settings,'LoopFieldName2') || strcmp(settings.LoopFieldName2,'NaN') || strcmp(settings.LoopFieldName2,'Coil_Cycle_Order')
        % Correlation plots
        xlimmin = 0; xlimmax = 1;
        ylimmin = 0; ylimmax = 1;
        Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
        Alphabet = {'i','ii','iii','iv'};
        fontsize = 10; Font = 'calibri'; Axis_fontsize = 11;
        
        if size(settings.LoopValues2,1) == 3
            figure('color','w','Units','normalized','Position',[0.03645833333333,0.0462962962963,0.3625,0.991666666666667])
        else
            figure('color','w','Units','normalized','Position',[0.085416666666667,0.075925925925926,0.384375,0.673148148148148])
        end
        %cmap = jet(8);
        cmap = [255,0,0;0,176,240]./255;
        tiledlayout('flow','tilespacing','compact','padding','none');
        n = 1; text_inc = 0.07;
        for Additional_Loop_Counter = 1:size(settings.LoopValues,1)
            for Additional_Loop_Counter2 = 1:size(settings.LoopValues2,1)
                if strcmpi(settings.LoopFieldName2,'Coil_Cycle_Order') && all(settings.LoopValues2(Additional_Loop_Counter2,:) == 1:8)
                    titlepart2 = 'Sequential ';
                elseif strcmpi(settings.LoopFieldName2,'Coil_Cycle_Order') && all(settings.LoopValues2(Additional_Loop_Counter2,:) == [1,4,7,2,5,8,3,6])
                    titlepart2 = 'Distributed ';
                else
                    titlepart2 = ' ';
                end
                fulltitle = [titlepart1{Additional_Loop_Counter},titlepart2];
                
                nexttile; textx = xlimmin + text_inc; texty = ylimmax - text_inc;
                for Tx_n = 1:8
                    plot(abs(Rel_Long_Maps(:,Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2)),abs(GT_Rel_Long_Maps(:,Tx_n,Additional_Loop_Counter,Additional_Loop_Counter2)),'color',cmap(Additional_Loop_Counter2,:),'Marker','.','Linestyle','none','Markersize',2,'HandleVisibility','off'); hold on; grid on
                end
                plot([xlimmin,xlimmax],[ylimmin,ylimmax],'k','HandleVisibility','off')
                
                Long_Maps = reshape(abs(Rel_Long_Maps(:,:,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,:,:)),[],size(settings.LoopValues,1),size(settings.LoopValues2,1));
                Long_GT_Maps = reshape(abs(GT_Rel_Long_Maps),[],size(settings.LoopValues,1),size(settings.LoopValues2,1));
                
                P = polyfit(Long_Maps(:,Additional_Loop_Counter,Additional_Loop_Counter2),Long_GT_Maps(:,Additional_Loop_Counter,Additional_Loop_Counter2),1);
                text(textx,texty,['y = ',num2str(P(1),'%4.2f'),'x ',num2strpn(P(2),'%4.2f')],'Fontsize',fontsize); texty = texty - text_inc;
                [~,~,~,~,stats] = regress(Long_Maps(:,Additional_Loop_Counter,Additional_Loop_Counter2),Long_GT_Maps(:,Additional_Loop_Counter,Additional_Loop_Counter2));
                text(textx,texty,['r^2 = ',num2str(stats(1),'%4.2f')],'Fontsize',fontsize);  texty = texty - text_inc;
                text(textx,texty,['n = ',num2str(size(Long_Maps,1))],'Fontsize',fontsize);  texty = texty - text_inc;
                text(textx,texty,['RMSE = ',num2str(mean(NRMSE_Values(:,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2),1),'%4.2f')],'Fontsize',fontsize);
                
                title(strjoin([Alphabet(n),') ',fulltitle]),'FontName',Font,'Fontsize',14)
                xlabel('Ground-Truth [a.u.]','FontName',Font,'FontSize',Axis_fontsize,'FontWeight','normal')
                ylabel('Relative B_1^+ [a.u.]','FontName',Font,'FontSize',Axis_fontsize,'FontWeight','normal')
                box on
                axis square
                yticks([0:0.2:1]);
                xticks([0:0.2:1]);
                xlim([xlimmin xlimmax])
                ylim([ylimmin ylimmax])
                
                % Legend
                %if Additional_Loop_Counter == 1 && Additional_Loop_Counter2 == 1
                %for Tx_n = 1:8
                %    plot(-10,-10,'color',cmap(Tx_n,:),'Marker','.','LineStyle','none','MarkerSize',20)
                %end
                %leg = legend('1','2','3','4','5','6','7','8','Location','southeast');
                %title(leg,'Tx #');
                %end
                n = n + 1;
            end
        end
    end
    
    
    if strcmp(settings.LoopFieldName2,'nom_FA')
        figure('color','w','Units','Normalized','Position',[0.490625,0.188888888888889,0.242708333333334,0.376851851851852])
        cmap = [0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.8500 0.3250 0.0980]; Markers = {'.','.','.'};
        for Additional_Loop_Counter = 1:size(settings.LoopValues,1)
            for Additional_Loop_Counter2 = 1:size(settings.LoopValues2,1)
                for Tx_n = 1:8
                    plot(settings.LoopValues2(Additional_Loop_Counter2).*180/pi,100.*NRMSE_Values(Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2),'color',cmap(Additional_Loop_Counter,:),'Marker',Markers{Additional_Loop_Counter},'MarkerSize',10,'Linewidth',2,'HandleVisibility','off'); hold on
                end
            end
        end
        xlim([0 settings.LoopValues2(end).*180/pi]); ylim([0 100]);
        grid on
        axis square
        ylabel('NRMSE, [%]')
        xlabel(['Nominal Flip Angle, [',char(176),']']);
        
        if strcmpi(settings.UseSyntheticData,'Duke')
            for Additional_Loop_Counter = 1:size(settings.LoopValues,1)
                for Additional_Loop_Counter2 = 1:size(settings.LoopValues2,1)
                    for Tx_n = 1:8
                        plot(settings.LoopValues2(Additional_Loop_Counter2).*180/pi,100.*Heart_NRMSE_Values(Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2),'color',cmap(Additional_Loop_Counter,:),'Marker','o','MarkerSize',6,'Linewidth',0.5,'HandleVisibility','off'); hold on
                    end
                end
            end
            p3 = plot(-10,-10,'color','k','Marker','.','LineStyle','none','MarkerSize',20,'Linewidth',1);
            p4 = plot(-10,-10,'color','k','Marker','o','LineStyle','none','MarkerSize',6,'Linewidth',0.5);
        end
        hold off
        text1 = ['\color[rgb]{',num2str(cmap(1,1)),' ',num2str(cmap(1,2)),' ',num2str(cmap(1,3)),'}',titlepart1{1}];
        text2 = ['\color[rgb]{',num2str(cmap(2,1)),' ',num2str(cmap(2,2)),' ',num2str(cmap(2,3)),'}',titlepart1{2}];
        annotation('textbox',[0.184358176612133,0.775593773133435,0.285407731487003,0.133906636366973],'String',[text1,newline,text2],'FitBoxToText','on','BackgroundColor','w','HorizontalAlignment','center');
        if strcmpi(settings.UseSyntheticData,'Duke')
            ah1 = axes('position',get(gca,'position'),'visible','off');
            leg2 = legend(ah1, [p3 p4], 'Body','Heart', 'Location',[0.708965485462648,0.788185910842733,0.144670960405183,0.112100739443917]); leg2.ItemTokenSize(1) = 10;
        end
    end
    
    if strcmp(settings.LoopFieldName2,'T1s') && settings.Global_T1 == 1
        figure('color','w','Units','Normalized','Position',[0.490625,0.188888888888889,0.242708333333334,0.376851851851852])
        cmap = [0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.8500 0.3250 0.0980]; Markers = {'.','.','.'};
        T1_n = 1;
        for Additional_Loop_Counter = 1:size(settings.LoopValues,1)
            for Additional_Loop_Counter2 = 1:size(settings.LoopValues2,1)
                for Tx_n = 1:8
                    plot(settings.LoopValues2(Additional_Loop_Counter2),100.*NRMSE_Values(Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2),'color',cmap(Additional_Loop_Counter,:),'Marker',Markers{Additional_Loop_Counter},'MarkerSize',10,'Linewidth',2,'HandleVisibility','off'); hold on
                end
            end
        end
        xlim([0 settings.LoopValues2(end)]); ylim([0 100]);
        grid on
        axis square
        ylabel('NRMSE, [%]')
        xlabel('T_1, [s]'); 
        
        if strcmpi(settings.UseSyntheticData,'Duke')
            for Additional_Loop_Counter = 1:size(settings.LoopValues,1)
                for Additional_Loop_Counter2 = 1:size(settings.LoopValues2,1)
                    for Tx_n = 1:8
                        plot(settings.LoopValues2(Additional_Loop_Counter2),100.*Heart_NRMSE_Values(Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2),'color',cmap(Additional_Loop_Counter,:),'Marker','o','MarkerSize',6,'Linewidth',0.5,'HandleVisibility','off'); hold on
                    end
                end
            end
            p3 = plot(-10,-10,'color','k','Marker','.','LineStyle','none','MarkerSize',20,'Linewidth',1);
            p4 = plot(-10,-10,'color','k','Marker','o','LineStyle','none','MarkerSize',6,'Linewidth',0.5);
        end
        hold off
        text1 = ['\color[rgb]{',num2str(cmap(1,1)),' ',num2str(cmap(1,2)),' ',num2str(cmap(1,3)),'}',titlepart1{1}];
        text2 = ['\color[rgb]{',num2str(cmap(2,1)),' ',num2str(cmap(2,2)),' ',num2str(cmap(2,3)),'}',titlepart1{2}];
        annotation('textbox',[0.197770193779515,0.775593773133434,0.285407731487003,0.133906636366973],'String',[text1,newline,text2],'FitBoxToText','on','BackgroundColor','w','HorizontalAlignment','center');
        if strcmpi(settings.UseSyntheticData,'Duke')
            ah1 = axes('position',get(gca,'position'),'visible','off');
            leg2 = legend(ah1, [p3 p4], 'Body','Heart', 'Location',[0.708965485462648,0.788185910842733,0.144670960405183,0.112100739443917]); leg2.ItemTokenSize(1) = 10;
        end
    end
    
    if strcmp(settings.LoopFieldName2,'Matrix_Size')
        figure('color','w','Units','Normalized','Position',[0.490625,0.188888888888889,0.242708333333334,0.376851851851852])
        %cmap = jet(8); Markers = {'x','o'};
        cmap = [0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.8500 0.3250 0.0980]; Markers = {'.','.','.'};
        for Additional_Loop_Counter = 1:size(settings.LoopValues,1)
            for Additional_Loop_Counter2 = 1:size(settings.LoopValues2,1)
                for Tx_n = 1:8
                    plot(settings.LoopValues2(Additional_Loop_Counter2),100.*NRMSE_Values(Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2),'color',cmap(Additional_Loop_Counter2,:),'Marker',Markers{Additional_Loop_Counter},'MarkerSize',10,'Linewidth',2,'HandleVisibility','off'); hold on
                end
            end
        end
        %       plot(-10,-10,'color','k','Marker',Markers{1},'LineStyle','none','MarkerSize',10,'Linewidth',2)
        %       plot(-10,-10,'color','k','Marker',Markers{2},'LineStyle','none','MarkerSize',10,'Linewidth',2)
        plot(-10,-10,'color',cmap(1,:),'Marker',Markers{1},'LineStyle','none','MarkerSize',20,'Linewidth',2)
        plot(-10,-10,'color',cmap(2,:),'Marker',Markers{2},'LineStyle','none','MarkerSize',20,'Linewidth',2)
        %plot(-10,-10,'color',cmap(3,:),'Marker',Markers{3},'LineStyle','none','MarkerSize',20,'Linewidth',2)
        legend(titlepart1{1},titlepart1{2},'Location','NorthWest');
        xlim([0 72]); ylim([0 100]);
        grid on
        axis square
        ylabel('NRMSE, [%]')
        xlabel('N_{PE}, [s]');
    end
    
    if strcmp(settings.LoopFieldName2,'IT_TR')
        figure('color','w','Units','Normalized','Position',[0.490625,0.188888888888889,0.242708333333334,0.376851851851852])
        cmap = [0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.8500 0.3250 0.0980]; Markers = {'.','.','.'};
        for Additional_Loop_Counter = 1:size(settings.LoopValues,1)
            for Additional_Loop_Counter2 = 1:size(settings.LoopValues2,1)
                for Tx_n = 1:8
                    plot(1e3.*settings.LoopValues2(Additional_Loop_Counter2),100.*NRMSE_Values(Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2),'color',cmap(Additional_Loop_Counter,:),'Marker',Markers{Additional_Loop_Counter},'MarkerSize',10,'Linewidth',2,'HandleVisibility','off'); hold on
                end
            end
        end
        xlim([0 1e3.*settings.LoopValues2(end)]); ylim([0 100]);
        grid on
        axis square
        ylabel('NRMSE, [%]')
        xlabel('TR, [ms]');
        
        if strcmpi(settings.UseSyntheticData,'Duke')
            for Additional_Loop_Counter = 1:size(settings.LoopValues,1)
                for Additional_Loop_Counter2 = 1:size(settings.LoopValues2,1)
                    for Tx_n = 1:8
                        plot(1e3.*settings.LoopValues2(Additional_Loop_Counter2),100.*Heart_NRMSE_Values(Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2),'color',cmap(Additional_Loop_Counter,:),'Marker','o','MarkerSize',6,'Linewidth',0.5,'HandleVisibility','off'); hold on
                    end
                end
            end
            p3 = plot(-10,-10,'color','k','Marker','.','LineStyle','none','MarkerSize',20,'Linewidth',1);
            p4 = plot(-10,-10,'color','k','Marker','o','LineStyle','none','MarkerSize',6,'Linewidth',0.5);
        end
        hold off
        text1 = ['\color[rgb]{',num2str(cmap(1,1)),' ',num2str(cmap(1,2)),' ',num2str(cmap(1,3)),'}',titlepart1{1}];
        text2 = ['\color[rgb]{',num2str(cmap(2,1)),' ',num2str(cmap(2,2)),' ',num2str(cmap(2,3)),'}',titlepart1{2}];
        annotation('textbox',[0.184358176612133,0.775593773133435,0.285407731487003,0.133906636366973],'String',[text1,newline,text2],'FitBoxToText','on','BackgroundColor','w','HorizontalAlignment','center');
        if strcmpi(settings.UseSyntheticData,'Duke')
            ah1 = axes('position',get(gca,'position'),'visible','off');
            leg2 = legend(ah1, [p3 p4], 'Body','Heart', 'Location',[0.708965485462648,0.788185910842733,0.144670960405183,0.112100739443917]); leg2.ItemTokenSize(1) = 10;
        end
    end
    
    if strcmp(settings.LoopFieldName2,'Dummy_RF')
        figure('color','w','Units','Normalized','Position',[0.490625,0.188888888888889,0.242708333333334,0.376851851851852])
        %cmap = jet(8); Markers = {'x','o'};
        cmap = [0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.8500 0.3250 0.0980]; Markers = {'.','.','.'};
        for Additional_Loop_Counter = 1:size(settings.LoopValues,1)
            for Additional_Loop_Counter2 = 1:size(settings.LoopValues2,1)
                for Tx_n = 1:8
                    plot(settings.LoopValues2(Additional_Loop_Counter2),100.*NRMSE_Values(Tx_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Additional_Loop_Counter,Additional_Loop_Counter2),'color',cmap(Additional_Loop_Counter,:),'Marker',Markers{Additional_Loop_Counter},'MarkerSize',10,'Linewidth',2,'HandleVisibility','off'); hold on
                end
            end
        end
        %       plot(-10,-10,'color','k','Marker',Markers{1},'LineStyle','none','MarkerSize',10,'Linewidth',2)
        %       plot(-10,-10,'color','k','Marker',Markers{2},'LineStyle','none','MarkerSize',10,'Linewidth',2)
        plot(-10,-10,'color',cmap(1,:),'Marker',Markers{1},'LineStyle','none','MarkerSize',20,'Linewidth',2)
        plot(-10,-10,'color',cmap(2,:),'Marker',Markers{2},'LineStyle','none','MarkerSize',20,'Linewidth',2)
        %plot(-10,-10,'color',cmap(3),'Marker',Markers{3},'LineStyle','none','MarkerSize',10,'Linewidth',2)
        legend(titlepart1{1},titlepart1{2},'Location','NorthWest');
        xlim([0 settings.LoopValues2(end)]); ylim([0 100]);
        grid on
        axis square
        ylabel('NRMSE, [%]')
        xlabel('# Dummy RF (Per Channel), [a.u.]');
        %
    end
    
    % Plot SNR
    if strcmpi(settings.LoopFieldName,'Coil_Cycle') && (strcmpi(settings.LoopValues(Additional_Loop_Counter),'SS') || strcmpi(settings.LoopValues(Additional_Loop_Counter),'CC'))
        %     for Additional_Loop_Counter2 = 1:size(settings.LoopValues2,1)
        %         Signal_Ratio(:,:,:,:,:,:,:,:,:,1,Additional_Loop_Counter2) = abs(Rel_Images(:,:,:,:,:,:,:,:,:,find(strcmpi(settings.LoopValues,'CC'),1),Additional_Loop_Counter2))./abs(Rel_Images(:,:,:,:,:,:,:,:,:,find(strcmpi(settings.LoopValues,'SW'),1),Additional_Loop_Counter2));
        %     end
        %     figure('color','w','Units','Normalized','Position',[0.490625,0.188888888888889,0.242708333333334,0.376851851851852])
        %     histogram(Signal_Ratio(:,:,:,:,:,:,:,:,:,1,1),-1.3:.01:1.3,'facealpha',.5,'edgecolor','none')
        %     hold on
        %     histogram(Signal_Ratio(:,:,:,:,:,:,:,:,:,1,2),-1.3:.01:1.3,'facealpha',.5,'edgecolor','none')
        %     box off
        %     axis tight
        %     legend boxoff
    end
end


end