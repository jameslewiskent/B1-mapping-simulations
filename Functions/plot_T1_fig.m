function plot_T1_fig(results,settings,plot_settings)
% T1 plot

if strcmpi(settings.Scheme,'AFI')
    Nominal_FA = settings.Dynamic_Range.*settings.nomFA*(180/pi);
else
    Nominal_FA = settings.Dynamic_Range.*settings.nomPP_FA*(180/pi);
end

[Dynamic_Range,Dynamic_Range_Values] = Calc_Dynamic_Range(results,settings,plot_settings);
if plot_settings.Dynamic_Range_Axis == 1
    Dynamic_Range_Value = Dynamic_Range_Values(1);
    Axis_Values = Nominal_FA./Dynamic_Range_Values(1); % Rescale axis based on dynamic range
else
    Dynamic_Range_Value = 1;
    Axis_Values = Nominal_FA;
end

if plot_settings.Plot_Difference == 1
    Subtract_Linear = Axis_Values;
else
    Subtract_Linear = 0;
    plot([0,Axis_Values(end)],[0,Axis_Values(end)],'color',[0.8 0.8 0.8],'handlevisibility','off'); hold on
end

if size(settings.T1s,2) == 6
    cmap = sixcolourmap;
    choose_T1n = 1:6;
else
    cmap = turbo(size(settings.T1s,2));
    choose_T1n = round(linspace(1,size(settings.T1s,2),6));% pick 6 T1's to put in legend
end
Styles = {'-','--','-.',':'};
Noise_n = 1;
B0_n = 1;
Flow_n = 1;
Diff_n = 1;
for T1_n = 1:size(settings.T1s,2)
    for Repeat_n = 1:size(settings.Repeats,2)
        plot(Axis_Values,squeeze(results.Measured_FA(1,1,1,:,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Repeat_n))/Dynamic_Range_Value - Subtract_Linear','color',cmap(T1_n,:),'LineStyle','none','marker','.','markersize',3,'handlevisibility','off')
        hold on
    end
    h = plot(Axis_Values,squeeze(mean(results.Measured_FA(1,1,1,:,B0_n,T1_n,Flow_n,Diff_n,Noise_n,:),10))/Dynamic_Range_Value - Subtract_Linear','color',cmap(T1_n,:),'linewidth',1,'LineStyle',Styles{B0_n},'handlevisibility','on');
    if all(T1_n ~= choose_T1n)
        h.HandleVisibility = 'off';
    end
end
grid on
axis square
lgd = legend(sprintfc('%g', settings.T1s(choose_T1n)),'Location','NorthWest','Orientation','vertical');
lgd.Title.String = 'T1 (s)';
title(settings.title_string,'Fontsize',16)



if plot_settings.Dynamic_Range_Axis == 1 && plot_settings.Plot_Difference == 0
    xlabel('Applied Dynamic Range, [a.u.]'); xticks(1:10)
    ylabel('Measured Dynamic Range, [a.u.]');
    xlim([0 10]); ylim([0 10]);
elseif plot_settings.Dynamic_Range_Axis == 1 && plot_settings.Plot_Difference == 1
    xlabel('Applied Dynamic Range, [a.u.]'); xticks(1:10)
    ylabel('[Measured - Applied] Dynamic Range, [a.u.]');
    xlim([0 10]); ylim([-1 1]);
elseif plot_settings.Dynamic_Range_Axis == 0 && plot_settings.Plot_Difference == 0
    xlabel(['Nominal FA, [',char(176),']']); xticks(-20:20:220)
    ylabel(['Measured FA, [',char(176),']']);
    xlim([0 200]); ylim([0 200]);
    % Add dynamic range bar if there are sufficient repeats
    if settings.Repeats >= 20
        text((Dynamic_Range_Values(2)+Dynamic_Range_Values(1))/2,10,['DR: ',num2str(Dynamic_Range,'%.1f')],'horizontalalignment','center','verticalalignment','middle','color','k')
        plot(linspace(Dynamic_Range_Values(1),(Dynamic_Range_Values(2)+Dynamic_Range_Values(1))/2 -15,10),10*ones(1,10),'color','k','Linewidth',1,'handlevisibility','off')
        plot(linspace((Dynamic_Range_Values(2)+Dynamic_Range_Values(1))/2 +15,Dynamic_Range_Values(2),10),10*ones(1,10),'color','k','Linewidth',1,'handlevisibility','off')
        plot(Dynamic_Range_Values(1)*ones(1,10),linspace(8,12,10),'color','k','Linewidth',1,'handlevisibility','off')
        plot(Dynamic_Range_Values(2)*ones(1,10),linspace(8,12,10),'color','k','Linewidth',1,'handlevisibility','off')
    end
elseif plot_settings.Dynamic_Range_Axis == 0 && plot_settings.Plot_Difference == 1
    xlabel(['Nominal FA, [',char(176),']']); xticks(-20:20:220)
    ylabel(['[Measured - Nominal] FA, [',char(176),']']);
    xlim([0 200]); ylim([-100 100]);
    % Add dynamic range bar if there are sufficient repeats
    if settings.Repeats >= 20
        text((Dynamic_Range_Values(2)+Dynamic_Range_Values(1))/2,10,['DR: ',num2str(Dynamic_Range,'%.1f')],'horizontalalignment','center','verticalalignment','middle','color','k')
        plot(linspace(Dynamic_Range_Values(1),(Dynamic_Range_Values(2)+Dynamic_Range_Values(1))/2 -15,10),10*ones(1,10),'color','k','Linewidth',1,'handlevisibility','off')
        plot(linspace((Dynamic_Range_Values(2)+Dynamic_Range_Values(1))/2 +15,Dynamic_Range_Values(2),10),10*ones(1,10),'color','k','Linewidth',1,'handlevisibility','off')
        plot(Dynamic_Range_Values(1)*ones(1,10),linspace(8,12,10),'color','k','Linewidth',1,'handlevisibility','off')
        plot(Dynamic_Range_Values(2)*ones(1,10),linspace(8,12,10),'color','k','Linewidth',1,'handlevisibility','off')
    end
end

end

