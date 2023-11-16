function [results,settings] = run_sequence_simulations(settings,results)
% Function handles iterating through scheme(s) to simulate and hands over to analysing
% data in analysis_function_core.
%
% James Kent. 2023. Using B. Hargreaves EPG Functions.

settings.Scan_Size(1) = round(settings.PE1_Resolution*settings.PE1_Partial_Fourier*settings.Matrix_Size(1));
settings.Scan_Size(2) = round(settings.PE2_Resolution*settings.PE2_Partial_Fourier*settings.Matrix_Size(2));
settings = Calc_Slice_Shifts(settings); % Calculate Slice Shifts for MS sequences
settings.Mag_Track_Flags = zeros(1,size(settings.Mag_Track_FAValues,2));
settings.filepath = fullfile('Data',lower(settings.Scheme));
settings.lookup_filename = [settings.Scheme,'_lookup_table.mat'];

if ~any(settings.B0_Range_Hz == 0)
    disp('We require to simulate B0 = 0 Hz for lookup table calculation. Adding 0 Hz to B0_Range_Hz array.')
    settings.B0_Range_Hz = [0,settings.B0_Range_Hz];
end

if ~any(settings.Velocities == 0)
    disp('We require to simulate flow velocity  = 0 ms^{-1} for lookup table calculation. Adding 0 ms^{-1} to Velocities array.')
    settings.Velocities = [0,settings.Velocities];
    settings.Angles = [0,settings.Angles];
end

if ~any(settings.Diff_coeffs == 0)
    disp('We require to simulate diffusion coefficient = 0 m^{2}s^{-1} for lookup table calculation. Adding 0 m^{2}s^{-1} to Diff_coeffs array.')
    settings.Diff_coeffs = [0,settings.Diff_coeffs];
end

if ~any(isnan(settings.Noise))
    disp('We require to simulate zero noise for lookup table calculation. Adding no noise to Noise array.')
    settings.Noise = [NaN,settings.Noise];
end

if ~exist(settings.filepath,'dir')
mkdir(settings.filepath);
end
    
% Load sequence specific settings
settings = load_sequence_settings(settings);

% Calculate slice ordering
settings = Calc_Slice_Order(settings);

% Check if simulations already ran
[already_ran,filename] = check_if_already_ran(settings);

if already_ran % Don't re-simulate if input parameters unchanged
    disp('Input parameters unchanged - results not re-simulated.')
    load(fullfile(settings.filepath,filename),'results','settings')  
else % Simulate sequence
    disp('Starting simulations.')
    if strcmpi(settings.Scheme,'SatTFL')
        [simulation_results,settings] = epg_sattfl(settings);
    elseif strcmpi(settings.Scheme,'Sandwich')
        [simulation_results,settings] = epg_sandwich(settings);
    elseif strcmpi(settings.Scheme,'SA2RAGE')
        [simulation_results,settings] = epg_sa2rage(settings);
    elseif strcmpi(settings.Scheme,'AFI')
        [simulation_results,settings] = epg_afi(settings);
    elseif strcmpi(settings.Scheme,'DREAM')
        [simulation_results,settings] = epg_dream(settings);
    else
        error('ABORTED: Scheme not recognised, please input either ''SatTFL'', ''Sandwich'', ''DREAM'', ''AFI'', ''SA2RAGE'' OR ''ALL''.')
    end
    disp('Simulations successful!')
    
    % Process simulation results
    % Call for however modes are required
    disp('Starting analysis.')
    [results] = analysis_function_core(simulation_results, settings,results);
    disp('Analysis finished.')
    filename = ['Results_',char(datetime('now','TimeZone','local','Format','d-MMM-y-HH-mm')),'.mat'];
    save(fullfile(settings.filepath,filename),'results','settings')     
end

% Display some energy statistics
Pulses_string = 'Consisting of ';
if isfield(results,'N_imaging_RF') && results.N_imaging_RF ~= 0
Pulses_string =  [Pulses_string,num2str(results.N_imaging_RF),' imaging RF pulses'];
end
if isfield(results,'N_imaging_RF1') && results.N_imaging_RF1 ~= 0
Pulses_string =  [Pulses_string,num2str(results.N_imaging_RF1),' imaging RF1 pulses'];
end
if isfield(results,'N_imaging_RF2') && results.N_imaging_RF2 ~= 0
Pulses_string =  [Pulses_string,', ',num2str(results.N_imaging_RF2),' imaging RF2 pulses'];
end
if isfield(results,'N_prep_RF') && results.N_prep_RF ~= 0
Pulses_string = [Pulses_string,', ',num2str(results.N_prep_RF),' preparation RF pulses'];
end
disp([upper(settings.Scheme),' acquisition time = ',num2str(results.Cumulative_Time),' s. ',Pulses_string,'.'])
if ~isempty(find(settings.Dynamic_Range == 1, 1)) % e.g. nominal flip angle
    disp([upper(settings.Scheme),' sequence total energy = ',num2str(results.Total_Energy(settings.Dynamic_Range == 1)),' Joules (per channel) at nominal flip angle.'])
    disp([upper(settings.Scheme),' sequence average 10 second power = ',num2str(results.Average_10s_Power(settings.Dynamic_Range == 1)),' Watts (per channel) at nominal flip angle.'])
end

end


