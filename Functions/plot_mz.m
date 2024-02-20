function  plot_mz(results,settings)

cmap = get(gca,'colororder');
cmap = [cmap;0.8,0.8,0.8];
for n = 1:size(results.Mag_Track,2)
%     if n ~= 6
%         color_n = 8;
%     else
         color_n = n;
%     end
    
plot(results.Mag_Track{n}(2,:),results.Mag_Track{n}(1,:),'linewidth',2,'color',cmap(color_n,:)); hold on
%[x_ss] = Calc_SS(results.Mag_Track{n}(1,:));
%xline(results.Mag_Track{n}(2,x_ss))
end
xlabel('Time, [s]');
ylabel('Longitudinal Magnetisation, M_z')
ylim([-1 1]); yticks([-1 -0.5 0 0.5 1])
xlim([0 results.Mag_Track{n}(2,end)])

if length(settings.Mag_Track_FAValues) <= 3
    nCols = 1;
else
    nCols = 2;
end

if settings.UseSyntheticData == 0
lgd = legend(sprintfc('%g', settings.Mag_Track_FAValues),'Location','South','Orientation','vertical','NumColumns',nCols);
lgd.Title.String = ['Nominal FA, [',char(176),']'];
end
end

