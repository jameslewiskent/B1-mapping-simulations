function plot_B0_fig(results,settings)
% B0  plot
if strcmp(settings.Scheme,'AFI')
    Nominal_FA = settings.Dynamic_Range.*settings.nomFA*(180/pi);
else
    Nominal_FA = settings.Dynamic_Range.*settings.nomPP_FA*(180/pi);
end

[~,Dynamic_Range_Values] = Calc_Dynamic_Range(results,settings);
    
Dynamic_Range_Axis = Nominal_FA./Dynamic_Range_Values(1); % Rescale axis based on dynamic range
    

    cmap = sixcolourmap; 
    Styles = {'-','--','-.',':'};
Noise_n = 1;
ah = axes;
plot([0,200],[0,200],'color',[0.8 0.8 0.8],'handlevisibility','off'); hold on
T1_n = 3;
for B0_n = 1:size(settings.B0_Range_Hz,2)
    plot(Dynamic_Range_Axis,squeeze(mean(results.Measured_FA(1,1,1,:,B0_n,T1_n,Noise_n,:),8)),'color',cmap(T1_n,:),'linewidth',1.5,'LineStyle',Styles{B0_n})
    hold on
end
grid on
axis square
lgd = legend(ah,sprintfc('%g', settings.B0_Range_Hz/1e3),'Location','NorthWest','Orientation','vertical');
lgd.Title.String = '\Deltaf_0 (kHz)';
xlabel('Applied Dynamic Range, [a.u.]'); xticks(1:10)
ylabel(['Measurable Dynamic Range, [',char(176),']']);
xlim([0 10]); ylim([0 10]);
title(settings.title_string,'Fontsize',16)
end

