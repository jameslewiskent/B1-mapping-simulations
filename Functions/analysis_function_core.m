function [results] = analysis_function_core(simulation_results, settings,results)
% Function handles the actually analysis i.e. Reordering, zero-filling,
% generating synthetic noise, Fourier transform, calculating alpha.
disp('Starting analysis.')

[results.Max_Val_IT1,FWHM1] = Process_IT(settings,simulation_results.IT1,1);
results.FWHM = FWHM1;
if isfield(simulation_results,'IT2')
    [results.Max_Val_IT2,FWHM2] = Process_IT(settings,simulation_results.IT2,2);
    results.FWHM = cat(1,FWHM1,FWHM2);
else
    results.Max_Val_IT2 = NaN;
end
if settings.T1Corr == 1
    [results.Max_Val_IT3,FWHM3] = Process_IT(settings,simulation_results.IT3,3);
    results.FWHM = cat(1,FWHM1,FWHM2,FWHM3);
else
    results.Max_Val_IT3 = NaN;
end

% Calculate FA
[results.Measured_FA] = Calc_FA(settings,results.Max_Val_IT1,results.Max_Val_IT2,results.Max_Val_IT3);

if (strcmpi(settings.UseSyntheticData,'Duke') || strcmpi(settings.UseSyntheticData,'Phantom'))
    % Reshape long synthetic data to image (might be a single mode (e.g. CP) or multitransmit)
    Image_Maps = reshape(permute(results.Measured_FA,[4,2,3,1,5,6,7,8,9]),[size(settings.Tx_FA_map),1,size(results.Measured_FA,5:ndims(results.Measured_FA))]);
    
    % Perform unencoding
    if settings.Modes > 1 && ~strcmpi(settings.Enc_Scheme,'Indiv')
        Image_Maps = Pixelwise_Unencoding(settings,Image_Maps);
    end
    
    if strcmpi(settings.Scheme,'GRE')
        % Calculate relative maps
        results.Image_Maps = Image_Maps./sum(abs(Image_Maps),3);
    else
        results.Image_Maps = Image_Maps;
    end
end

% Pass simulation results
if isfield(simulation_results,'seq_TRs')
    results.seq_TRs = simulation_results.seq_TRs;
end
results.Total_Energy = simulation_results.Total_Energy;
results.Average_10s_Power = simulation_results.Average_10s_Power;
results.Mag_Track = simulation_results.Mag_Track;
results.Cumulative_Time = simulation_results.Cumulative_Time;
if isfield(simulation_results,'N_imaging_RF')
    results.N_imaging_RF = simulation_results.N_imaging_RF;
end
if isfield(simulation_results,'N_prep_RF')
    results.N_prep_RF = simulation_results.N_prep_RF;
end
if isfield(simulation_results,'N_imaging_RF1')
    results.N_imaging_RF1 = simulation_results.N_imaging_RF1;
end
if isfield(simulation_results,'N_imaging_RF2')
    results.N_imaging_RF2 = simulation_results.N_imaging_RF2;
end

disp('Analysis finished.')

end

