function  plot_mz(results,settings)

% cmap = get(gca,'colororder');
% cmap = [cmap;0.8,0.8,0.8];
cmap = jet(size(results.Mag_Track,2)); flag = 0;
for n = 1:size(results.Mag_Track,2)    
plot(results.Mag_Track{n}(2,:),results.Mag_Track{n}(1,:),'linewidth',2,'color',cmap(n,:)); hold on

% Change axis depending on lowest Mz
if any(results.Mag_Track{n}(1,:) <0)
   ylimmin = -1; 
   flag = 1;
elseif flag == 0
   ylimmin = 0;
end
end



xlabel('Time, [s]');
ylabel('Longitudinal Magnetisation, M_z')
ylim([ylimmin 1]); yticks(ylimmin:0.5:1)
xlim([0 results.Mag_Track{n}(2,end)])

if length(settings.Mag_Track_FAValues) <= 3
    nCols = 1;
else
    nCols = 2;
end

if ~(strcmpi(settings.UseSyntheticData,'Duke') || strcmpi(settings.UseSyntheticData,'Phantom'))
lgd = legend(sprintfc('%g', settings.Mag_Track_FAValues),'Location','SouthWest','Orientation','vertical','NumColumns',nCols);
lgd.Title.String = ['Nominal FA, [',char(176),']'];
end
end

