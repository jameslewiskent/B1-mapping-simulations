function plot_Diff_fig(results,settings,plot_settings)
% Diffusion  plot

if strcmpi(settings.Scheme,'AFI')
    Nominal_FA = settings.Dynamic_Range.*settings.nomFA*(180/pi);
else
    Nominal_FA = settings.Dynamic_Range.*settings.nomPP_FA*(180/pi);
end

if plot_settings.Dynamic_Range_Axis == 1
    [~,Dynamic_Range_Values] = Calc_Dynamic_Range(results,settings,plot_settings);
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
Noise_n = find(~isnan(settings.Noise)); % First non-NaN (non-zero noise)
if isempty(Noise_n)
    Noise_n = 1;
else
    Noise_n = Noise_n(1);
end
B0_n = 1;
T1_n = 3;
Flow_n = 1;
for Diff_n = 1:size(settings.Diff_coeffs,2)
    plot(Axis_Values,squeeze(mean(results.Measured_FA(1,1,1,:,B0_n,T1_n,Flow_n,Diff_n,Noise_n,:),10))/Dynamic_Range_Value - Subtract_Linear','color',cmap(T1_n,:),'linewidth',1.5,'LineStyle',Styles{Diff_n});    hold on
end
grid on
axis square
lgd = legend(sprintfc('%g', settings.Diff_coeffs),'Location','NorthWest','Orientation','vertical');
lgd.Title.String = 'Diffusion Coeff. (m^2s^{-1})';
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

