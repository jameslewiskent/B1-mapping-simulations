function plot_FWHM(results,settings,plot_settings)
B0_n = 1;
T1_n = 3;
Flow_n = 1;
Diff_n = 1;

figure('Name','FWHMs','color','w','Position',[488,342,776.2,420])
tiledlayout(1,2,'tilespacing','compact','padding','none');

for Scheme_n = 1:size(results,2)
    Noise_n = find(~isnan(settings(1,Scheme_n).Noise)); % First non-NaN (non-zero noise)
if isempty(Noise_n)
    Noise_n = 1;
else
    Noise_n = Noise_n(1);
end

    if strcmpi(settings(1,Scheme_n).Scheme,'AFI')
        Nominal_FA = settings(1,Scheme_n).Dynamic_Range.*settings(1,Scheme_n).nom_FA*(180/pi);
    else
        Nominal_FA = settings(1,Scheme_n).Dynamic_Range.*settings(1,Scheme_n).nomPP_FA*(180/pi);
    end
    if plot_settings.Dynamic_Range_Axis == 1
        Axis_Values(:,Scheme_n) = settings(1,Scheme_n).Dynamic_Range;
    else
        Axis_Values(:,Scheme_n) = Nominal_FA;
    end
    
    nexttile(1);
    if size(settings,2) > 1
    plot(Axis_Values(:,Scheme_n),squeeze(mean(abs(results(1,Scheme_n).FWHM(1,1,:,B0_n,T1_n,Flow_n,Diff_n,Noise_n,:)),9)),'linewidth',1.5); hold on % try mesh or imagesc
    else
    plot(Axis_Values(:,Scheme_n),squeeze(mean(abs(results(1,Scheme_n).FWHM(1,1,:,B0_n,:,Flow_n,Diff_n,Noise_n,:)),[5,9])),'r.','HandleVisibility','off'); hold on % try mesh or imagesc
    if any(results(1,Scheme_n).FWHM) ~= 0
    plot(Axis_Values(:,Scheme_n),squeeze(mean(abs(results(1,Scheme_n).FWHM(1,2,:,B0_n,:,Flow_n,Diff_n,Noise_n,:)),[5,9])),'b.','HandleVisibility','off');
    end
    end
    caxis([0 7])
    title('S_0')
    ylabel('PSF FWHM, [voxels]');
    ylim([1 4])
    if plot_settings.Dynamic_Range_Axis == 1
    xlim([0 max(Axis_Values,[],'all')])
    xlabel('Applied Dynamic Range, [a.u.]');
    else
    xlim([0 200])
    xlabel(['Nominal Presaturation Flip Angle, [',char(176),']']);
    end
    axis square
    grid on
    
    nexttile(2);
    if size(settings,2) > 1
    plot(Axis_Values(:,Scheme_n),squeeze(mean(abs(results(1,Scheme_n).FWHM(2,1,:,B0_n,T1_n,Flow_n,Diff_n,Noise_n,:)),9)),'linewidth',1.5); hold on
    else
    plot(Axis_Values(:,Scheme_n),squeeze(mean(abs(results(1,Scheme_n).FWHM(2,1,:,B0_n,:,Flow_n,Diff_n,Noise_n,:)),[5,9])),'r.','HandleVisibility','off'); hold on
    if any(results(1,Scheme_n).FWHM)
    plot(Axis_Values(:,Scheme_n),squeeze(mean(abs(results(1,Scheme_n).FWHM(2,2,:,B0_n,:,Flow_n,Diff_n,Noise_n,:)),[5,9])),'b.','HandleVisibility','off'); hold on
    end
    end
    caxis([0 7])
    title('S_1')
    ylabel('PSF FWHM, [voxels]');
    ylim([1 4])
    if plot_settings.Dynamic_Range_Axis == 1
    xlim([0 max(Axis_Values,[],'all')])
    xlabel('Applied Dynamic Range, [a.u.]');
    else
    xlim([0 200])
    xlabel(['Nominal Presaturation Flip Angle, [',char(176),']']);
    end
    axis square
    grid on
end
if size(settings,2) > 1
nexttile(1);
lgdstr = {'SatTFL','Sandwich','SA2RAGE','AFI','DREAM'};
else
    plot(-1,-1,'r.','markersize',10); hold on
    plot(-1,-1,'b.','markersize',10)
    if any(results(1,Scheme_n).FWHM)
        lgdstr = {'PE1','PE2'};
    else
        lgdstr = {'PE1'};
    end
end
legend(lgdstr,'Location','NorthEast')


end