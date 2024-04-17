function plot_SAR_fig(results,settings)
figure('color','w');
for Scheme_n = 1:length(settings)
    
    if strcmpi(settings(1,Scheme_n).Scheme,'AFI')
        Nominal_FA = settings(1,Scheme_n).Dynamic_Range.*settings(1,Scheme_n).nom_FA*(180/pi);
    else
        Nominal_FA = settings(1,Scheme_n).Dynamic_Range.*settings(1,Scheme_n).nomPP_FA*(180/pi);
    end
    plot(Nominal_FA, results(1,Scheme_n).Average_10s_Power,'Linewidth',2); hold on
    lgdstr(Scheme_n) = {settings(1,Scheme_n).Scheme};
end

grid on
axis square
lgd = legend(string(lgdstr),'Location','NorthWest','Orientation','vertical');
xlabel(['Nominal FA, [',char(176),']']); xticks([-20:20:220])
ylabel(['Average 10s Power, [Watts]']);
xlim([0 200]); ylim([0 20]);
title('Average 10s Power','Fontsize',16)
end

