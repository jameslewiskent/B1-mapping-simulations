function plot_Additional_Loop_fig(results,settings,plot_settings)
% Additional Loop Plots

% Plot lookup table for each TR
if strcmp(settings.LoopFieldName,'TR')
    % Load in various files from the loop
    TR_Value = 1; % TR value to compare lookup tables against
    
    for Additional_Loop_Counter = 1:length(settings.LoopValues)
        lookup_filename = [settings.Scheme,'_lookup_table_Loop',num2str(Additional_Loop_Counter),'.mat'];
        load([settings.filepath,filesep,lookup_filename],'x_query','fx_interp');
        fx_interps(Additional_Loop_Counter,:) = fx_interp;
        x_queries(Additional_Loop_Counter,:) = x_query;
    end
    
    cmap = turbo(length(settings.LoopValues));
    figure('color','w','Units','normalized','position',[0.2,0.3,0.575,0.486]); tiledlayout('flow','padding','none','tilespacing','compact')
    nexttile;
    for Additional_Loop_Counter = 1:length(settings.LoopValues)
        plot(x_queries(Additional_Loop_Counter,:).*(180/pi),fx_interps(Additional_Loop_Counter,:),'color',cmap(Additional_Loop_Counter,:)); hold on
    end
    xlim([0 180])
    ylim([-1 1])
    title('a) Lookup Tables')
    lgd = legend(strsplit(num2str(settings.LoopValues)),'numcolumns',4,'location','southwest');
    lgd.Title.String = 'TR (s)';
    lgd.ItemTokenSize(1) = 10;
    ylabel('Image Ratio, [a.u.]'); xlabel(['Nominal Flip Angle, [',char(176),']'])
    grid on
    
    nexttile;
    for Additional_Loop_Counter = 1:length(settings.LoopValues)
        % Apply lookup tables
        for index = 1:length(x_queries(Additional_Loop_Counter,:))
            [~,min_ind] = min(abs(fx_interps(Additional_Loop_Counter,index) - fx_interps(settings.LoopValues == TR_Value,:)));
            error(Additional_Loop_Counter,index) = x_queries(Additional_Loop_Counter,min_ind) - x_queries(Additional_Loop_Counter,index);
        end
        % Plot error
        plot(x_queries(Additional_Loop_Counter,:).*(180/pi),error(Additional_Loop_Counter,:).*(180/pi),'color',cmap(Additional_Loop_Counter,:)); hold on
    end
    xlim([0 180])
    ylim([-10 10])
    title('a) Error in FA')
    %lgd = legend(strsplit(num2str(settings.LoopValues)),'numcolumns',4,'location','southwest');
    %lgd.Title.String = 'TR (s)';
    ylabel(['Error in Flip Angle, [',char(176),']']); xlabel(['Nominal Flip Angle, [',char(176),']'])
    grid on
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

end