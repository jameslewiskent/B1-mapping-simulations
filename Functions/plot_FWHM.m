function plot_FWHM(results,settings,plot_settings)
T1_n = 3;
figure('Name','FWHMs','color','w')
tiles = tiledlayout(1,2,'tilespacing','compact','padding','none');

for Scheme_n = 1:size(results,2)
    
    if strcmpi(settings(1,Scheme_n).Scheme,'AFI')
        Nominal_FA = settings(1,Scheme_n).Dynamic_Range.*settings(1,Scheme_n).nomFA*(180/pi);
    else
        Nominal_FA = settings(1,Scheme_n).Dynamic_Range.*settings(1,Scheme_n).nomPP_FA*(180/pi);
    end
    if plot_settings.Dynamic_Range_Axis == 1
        [~,Dynamic_Range_Values] = Calc_Dynamic_Range(results(1,Scheme_n),settings(1,Scheme_n),plot_settings);
        Axis_Values(:,Scheme_n) = Nominal_FA./Dynamic_Range_Values(1); % Rescale axis based on dynamic range
    else
        Dynamic_Range_Values = 1;
        Axis_Values(:,Scheme_n) = Nominal_FA;
    end
    
    nexttile(1); plot(Axis_Values(:,Scheme_n),squeeze(mean(abs(results(1,Scheme_n).FWHM(1,:,1,T1_n,1,:)),6)),'linewidth',1.5); hold on % try mesh or imagesc
    caxis([0 7])
    title('Image 1')
    ylabel('PSF FWHM, [a.u.]');
    xlabel('Applied Dynamic Range, [a.u.]');
    ylim([1.1 2])
    if plot_settings.Dynamic_Range_Axis == 1
    xlim([0 10])
    else
    xlim([0 200])    
    end
    
    nexttile(2); plot(Axis_Values(:,Scheme_n),squeeze(mean(abs(results(1,Scheme_n).FWHM(2,:,1,T1_n,1,:)),6)),'linewidth',1.5); hold on
    caxis([0 7])
    title('Image 2')
    ylabel('PSF FWHM, [a.u.]');
    xlabel('Applied Dynamic Range, [a.u.]');
    ylim([1.1 2])
    if plot_settings.Dynamic_Range_Axis == 1
    xlim([0 10])
    else
    xlim([0 200])    
    end
end
nexttile(1);
legend({'SatTFL','Sandwich','SA2RAGE','AFI','DREAM'},'Location','Northwest')


end