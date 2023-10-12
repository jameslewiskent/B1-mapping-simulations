function plot_Flow_fig(results,settings,plot_settings)
% Flow  plot

if strcmpi(settings.Scheme,'AFI')
    Nominal_FA = settings.Dynamic_Range.*settings.nomFA*(180/pi);
else
    Nominal_FA = settings.Dynamic_Range.*settings.nomPP_FA*(180/pi);
end

[~,Dynamic_Range_Values] = Calc_Dynamic_Range(results,settings,plot_settings);
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

cmap = sixcolourmap;
Styles = {'-','--','-.',':'};
Noise_n = 1;
B0_n = 1;
T1_n = 3;
Diff_n = 1;
for Flow_n = 1:size(settings.Velocities,2)
    plot(Axis_Values,squeeze(mean(results.Measured_FA(1,1,1,:,B0_n,T1_n,Flow_n,Diff_n,Noise_n,:),10))/Dynamic_Range_Value - Subtract_Linear','color',cmap(T1_n,:),'linewidth',1.5,'LineStyle',Styles{Flow_n});    hold on
end
grid on
axis square
lgd = legend(sprintfc('%g', settings.Velocities),'Location','NorthWest','Orientation','vertical');
lgd.Title.String = 'Flow Velocity (ms^{-1})';
title(settings.title_string,'Fontsize',16);

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
elseif plot_settings.Dynamic_Range_Axis == 0 && plot_settings.Plot_Difference == 1
    xlabel(['Nominal FA, [',char(176),']']); xticks(-20:20:220)
    ylabel(['[Measured - Nominal] FA, [',char(176),']']);
    xlim([0 200]); ylim([-100 100]);
end
end
