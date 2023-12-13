function [Dynamic_Range,first_last_values] = Calc_Dynamic_Range(results,settings,plot_settings)
% Function to find workable dynamic range
% Dynamic range defined as linear range where the difference in the mean
% flip angle and standard deviation is less than settings.percent_under% of the nominal value

if strcmpi(settings.Scheme,'AFI')
    x_axis = settings.Dynamic_Range.*settings.nomFA*(180/pi);
else
    x_axis = settings.Dynamic_Range.*settings.nomPP_FA*(180/pi);
end
x_axis = permute(x_axis,[4 3 1 2]);

Noise_n = find(~isnan(settings.Noise)); % Don't use NaN noise level

% Mean Range
Mean_FA = mean(results.Measured_FA(:,:,:,:,find(settings.B0_Range_Hz == 0),:,find(settings.Velocities == 0),find(settings.Diff_co == 0),Noise_n,:),[6,10]); % On resonance, average T1s and repeats
Mean_full_range = find(abs(Mean_FA - x_axis) < (x_axis*plot_settings.Dyn_Range_pc)); % remove gradient and find values below percent_under


% SD Range
SD_FA = std(results.Measured_FA(:,:,:,:,find(settings.B0_Range_Hz == 0),:,find(settings.Velocities == 0),find(settings.Diff_co == 0),Noise_n,:),[],[6,10]); % On resonance, average T1s and repeats
SD_full_range = find(abs(SD_FA) < (x_axis*plot_settings.Dyn_Range_pc));


first_last_range(1) = max([Mean_full_range(1),SD_full_range(1)]); % Take maximum of ranges
first_last_range(2) = min([Mean_full_range(end),find(diff(SD_full_range) >1,1)]);

%figure();plot(abs(Mean_FA - x_axis)); hold on; plot((x_axis*settings.Dyn_Range_pc))
%figure();plot(abs(SD_FA)); hold on; plot((x_axis*settings.Dyn_Range_pc))
% cmap = sixcolourmap;
% figure('color','w')
% for T1_n = 1:size(results.Measured_FA,4)
%     plot(x_axis,squeeze(mean(results.Measured_FA(:,:,find(settings.B0_Range_Hz == 0),T1_n,:,:),6)),'color',cmap(T1_n,:)); hold on
%     plot(x_axis,squeeze(std(results.Measured_FA(:,:,find(settings.B0_Range_Hz == 0),T1_n,:,:),[],6)),'--','color',cmap(T1_n,:)); hold on
% end
% plot(x_axis(first_last_range(1)),Mean_FA(first_last_range(1)),'rx');
% plot(x_axis(first_last_range(2)),Mean_FA(first_last_range(2)),'rx');

first_last_values(1) = x_axis(first_last_range(1));
first_last_values(2) = x_axis(first_last_range(2));

Dynamic_Range = first_last_values(2)./first_last_values(1);
end

