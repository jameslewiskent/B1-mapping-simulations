function plot_Additional_Loop_fig(results,settings,plot_settings)
% Additional Loop Plots
Noise_n = find(~isnan(settings.Noise)); % First non-NaN (non-zero noise)
if isempty(Noise_n)
    Noise_n = 1;
else
    Noise_n = Noise_n(1);
end

% Plot lookup table for each TR
if strcmp(settings.LoopFieldName,'TR')
    % Load in various files from the loop
    T1_n = 1:length(settings.T1s);
    B0_n = 1;
    Flow_n = 1;
    Diff_n = 1;
    for Additional_Loop_Counter = 1:length(settings.LoopValues)
        savefilename = [settings.filename,'_Loop',num2str(Additional_Loop_Counter),'.mat'];
        load([settings.filepath,filesep,savefilename],'results');
        FAs(Additional_Loop_Counter,:,:,:) = squeeze(mean(results.Measured_FA(1,1,1,:,B0_n,T1_n,Flow_n,Diff_n,Noise_n,:),[10]));
        DR(Additional_Loop_Counter) = results.Dynamic_Range;
        DR_Values(Additional_Loop_Counter,:) = results.DR_Values;
    end
    DR_ratio = DR./settings.LoopValues';
    figure('color','w','Units','normalized','position',[0.2,0.3,0.575,0.486]);
    plot(settings.LoopValues,DR,'rx'); hold on
    xlim([0.5 2])
    ylim([1 10])
    %ylim([floor(min(DR_ratio)*10)/10 ceil(max(DR_ratio)*10)/10])
    xticks(settings.LoopValues(1):0.2:settings.LoopValues(end))
    ylabel('Dynamic Range Efficiency, [s^{-1}]'); xlabel(['TR, [s]'])
    grid on
    axis square
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

if strcmp(settings.LoopFieldName,'Ejection_Fraction')
    Nominal_FA = settings.Dynamic_Range.*settings.nomPP_FA*(180/pi);
    lower_bound = Nominal_FA(35);
    upper_bound = Nominal_FA(160);
    T1_n = 1:length(settings.T1s);
    B0_n = 1;
    Flow_n = 1;
    Diff_n = 1;
    
    for Additional_Loop_Counter = 1:length(settings.LoopValues)
        filename = [settings.filename,'_Loop',num2str(Additional_Loop_Counter),'.mat'];
        load([settings.filepath,filesep,filename],'results');
        FAs(Additional_Loop_Counter,:,:,:) = squeeze(mean(results.Measured_FA(1,1,1,:,B0_n,T1_n,Flow_n,Diff_n,Noise_n,:),[10]));
    end
    
    Diff_FAs = (FAs(:,:,:) - Nominal_FA(1,:))./Nominal_FA(1,:); % Difference to EF = 0
    
    Mean_Diff_FAs = 100.*squeeze(mean(Diff_FAs,3));
    figure('color','w'); tiledlayout(1,2,'padding','none','tilespacing','compact');
    nexttile;
    imagesc(Nominal_FA,settings.LoopValues.*100,Mean_Diff_FAs,[-10 10]);
    set(gca,'YDir','normal')
    ylabel('Ejection Fraction, [%]');
    xlabel(['Nominal FA, [',char(176),']']);
    xlim([0 200]);
    axis square
    colormap(bluewhitered)
    cb = colorbar; cb.Label.String = 'Percent Error, [%]';
    xline(lower_bound,'linewidth',2,'linestyle','--');
    xline(upper_bound,'linewidth',2,'linestyle','--');
    title('a) Mean T1');
    
    Mean_Diff_FAs = 100.*squeeze(mean(Diff_FAs(:,lower_bound:upper_bound,:),2));
    nexttile();
    imagesc(settings.T1s,settings.LoopValues.*100,Mean_Diff_FAs,[-10 10]);
    set(gca,'YDir','normal')
    ylabel('Ejection Fraction, [%]');
    xlabel('T_1, [s]');
    axis square
    colormap(bluewhitered)
    cb = colorbar; cb.Label.String = 'Percent Error, [%]';
    title('b) Mean DR');
end

if strcmp(settings.LoopFieldName,'Coil_Cycle')
    T1_n = 1;
    B0_n = 1;
    Flow_n = 1;
    Diff_n = 1;
     for Additional_Loop_Counter = 1:length(settings.LoopValues)
        filename = [settings.filename,'_Loop',num2str(Additional_Loop_Counter),'.mat'];
        load([settings.filepath,filesep,filename],'results');
        Rel_Image_Maps(:,:,:,Additional_Loop_Counter) = squeeze(mean(results.Rel_Image_Maps(:,:,:,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,:),10));
     end
    
    % Calculate ground truth relative maps
    GT_Rel_Image_Maps = settings.IT_Tx_FA_map./sum(abs(settings.IT_Tx_FA_map),3);
    
    for Tx_n = 1:8
    % Remove NaNs
    Masked_GT_Rel_Image_Maps = settings.Synthetic_Mask.*abs(GT_Rel_Image_Maps(:,:,Tx_n));
    Masked_GT_Rel_Image_Maps(Masked_GT_Rel_Image_Maps == 0) = NaN;
    GT_Rel_Long_Maps = Masked_GT_Rel_Image_Maps(~isnan(Masked_GT_Rel_Image_Maps));
    
    Masked_NCC_Rel_Image_Maps = settings.Synthetic_Mask.*abs(Rel_Image_Maps(:,:,Tx_n,find(settings.LoopValues == 0)));
    NCC_Rel_Long_Maps = Masked_NCC_Rel_Image_Maps(~isnan(Masked_NCC_Rel_Image_Maps));

    Masked_CC_Rel_Image_Maps = settings.Synthetic_Mask.*abs(Rel_Image_Maps(:,:,Tx_n,find(settings.LoopValues == 1)));
    CC_Rel_Long_Maps = Masked_CC_Rel_Image_Maps(~isnan(Masked_CC_Rel_Image_Maps));
    
    NCC_nrmse(Tx_n) = rmse(NCC_Rel_Long_Maps,GT_Rel_Long_Maps ,'norm');
    CC_nrmse(Tx_n) = rmse(CC_Rel_Long_Maps,GT_Rel_Long_Maps ,'norm');
    end
    
    NCC_nrmse
    CC_nrmse
    
    figure('color','w')
    tiledlayout('flow','tilespacing','compact','padding','none'); nexttile;
    imagesc(imtile(abs(Rel_Image_Maps(:,:,:,find(settings.LoopValues == 1))),'GridSize', [2 4]),[0 1])
    axis image off
    nexttile;
    imagesc(imtile(abs(Rel_Image_Maps(:,:,:,find(settings.LoopValues == 0))),'GridSize', [2 4]),[0 1])
    axis image off
    
    % RMSE inside heart
    for Tx_n = 1:8
    % Remove NaNs
    Masked_GT_Rel_Image_Maps = settings.Whole_Heart_Mask.*abs(GT_Rel_Image_Maps(:,:,Tx_n));
    Masked_GT_Rel_Image_Maps(Masked_GT_Rel_Image_Maps == 0) = NaN;
    GT_Rel_Long_Maps = Masked_GT_Rel_Image_Maps(~isnan(Masked_GT_Rel_Image_Maps));
    
    Masked_NCC_Rel_Image_Maps = settings.Whole_Heart_Mask.*abs(Rel_Image_Maps(:,:,Tx_n,find(settings.LoopValues == 0)));
    Masked_NCC_Rel_Image_Maps(Masked_NCC_Rel_Image_Maps == 0) = NaN;
    NCC_Rel_Long_Maps = Masked_NCC_Rel_Image_Maps(~isnan(Masked_NCC_Rel_Image_Maps));

    Masked_CC_Rel_Image_Maps = settings.Whole_Heart_Mask.*abs(Rel_Image_Maps(:,:,Tx_n,find(settings.LoopValues == 1)));
    Masked_CC_Rel_Image_Maps(Masked_CC_Rel_Image_Maps == 0) = NaN;
    CC_Rel_Long_Maps = Masked_CC_Rel_Image_Maps(~isnan(Masked_CC_Rel_Image_Maps));
    
    NCC_nrmse(Tx_n) = rmse(NCC_Rel_Long_Maps,GT_Rel_Long_Maps ,'norm');
    CC_nrmse(Tx_n) = rmse(CC_Rel_Long_Maps,GT_Rel_Long_Maps ,'norm');
    end
    
    mean(NCC_nrmse)
    mean(CC_nrmse)
end


end