function [results] = analysis_function_core(simulation_results, settings,results)
% Function handles the actually analysis i.e. Reordering, zero-filling,
% generating synthetic noise, Fourier transform, calculating alpha.

[results.Max_Val_IT1,FWHM1] = Process_IT(settings,simulation_results.IT1,1);
if isfield(simulation_results,'IT2')
    [results.Max_Val_IT2,FWHM2] = Process_IT(settings,simulation_results.IT2,2);
else
    results.Max_Val_IT2 = NaN;
    FWHM2 = NaN(size(FWHM1));
end
if settings.T1Corr == 1
    [results.Max_Val_IT3,FWHM3] = Process_IT(settings,simulation_results.IT3,3);
else
    results.Max_Val_IT3 = NaN;
    FWHM3 = NaN(size(FWHM1));
end
results.FWHM = cat(1,FWHM1,FWHM2,FWHM3);

% Calculate FA
[results.Measured_FA] = Calc_FA(settings,results.Max_Val_IT1,results.Max_Val_IT2,results.Max_Val_IT3);

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
end

