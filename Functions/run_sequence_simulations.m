function [results,settings] = run_sequence_simulations(settings,results,plot_settings)
% Function handles iterating through scheme(s) to simulate and hands over to analysing
% data in analysis_function_core.
%
% James Kent. 2023. Using B. Hargreaves EPG Functions.

settings.Scan_Size(1) = round(settings.PE1_Resolution*settings.PE1_Partial_Fourier*settings.Matrix_Size(1));
settings.Scan_Size(2) = round(settings.PE2_Resolution*settings.PE2_Partial_Fourier*settings.Matrix_Size(2));
settings = Calc_Slice_Shifts(settings); % Calculate Slice Shifts for MS sequences

% Load sequence specific settings
settings = load_sequence_settings(settings);


if settings.UseSyntheticData == 1
[settings.Dynamic_Range,settings.Tx_FA_map,settings.Enc_Mat] = Calc_Tx(max(settings.RF_Pulse),settings); % Now pulse voltages are known, we can calculate the associated transmit field
disp(['Multi-transmit mode mapping is active. Simulating B1 Mapping of ', num2str(size(settings.Tx_FA_map,3)),' Transmit Mode Configuration.']);
end

settings.Mag_Track_Flags = zeros(1,size(settings.Mag_Track_FAValues,2));
settings.filepath = fullfile('Data',lower(settings.Scheme));
settings.lookup_filename = [settings.Scheme,'_lookup_table.mat'];

if settings.UseSyntheticData == 1 && settings.Modes < 8  && strcmp(settings.Enc_Scheme,'Indiv')
    disp('Multitransmit mapping is on, and encoding scheme is individual but the number of modes is less than the number of channels (8). Setting modes to 8.')
    settings.Modes = 8;
elseif settings.UseSyntheticData == 1 && settings.Modes > 8  && strcmp(settings.Enc_Scheme,'Indiv')
    disp('Multitransmit mapping is on, and encoding scheme is individual but the number of modes is more than the number of channels (8). Setting modes to 8.')
    settings.Modes = 8;
end

if settings.UseSyntheticData == 1
   settings.Mag_Track_FAValues = 1; 
end

if settings.UseSyntheticData == 1 && settings.Modes == 1 && ~strcmp(settings.Enc_Scheme,'CP')
    disp('For Modes = 1 encoding scheme is fixed as CP. Setting settings.Enc_Scheme = "CP"') 
    settings.Enc_Scheme = 'CP';     
end

if settings.UseSyntheticData == 1 && settings.Repeats > 1
    disp('Setting repeats to 1.')
   settings.Repeats = 1;
end

if ~any(settings.B0_Range_Hz == 0)
    disp('We require to simulate B0 = 0 Hz for lookup table calculation. Adding 0 Hz to B0_Range_Hz array.')
    settings.B0_Range_Hz = [0,settings.B0_Range_Hz];
end

if ~any(settings.Velocities == 0)
    disp('We require to simulate flow velocity  = 0 ms^{-1} for lookup table calculation. Adding 0 ms^{-1} to Velocities array.')
    settings.Velocities = [0,settings.Velocities];
    settings.Angles = [0,settings.Angles];
end

if length(settings.Velocities) ~= length(settings.Angles)
    disp('Velocities and Angles settings must have same number of entries')
   settings.Angles = repmat(settings.Angles(1),1,length(settings.Velocities));
end

if ~any(settings.Diff_coeffs == 0)
    disp('We require to simulate diffusion coefficient = 0 m^{2}s^{-1} for lookup table calculation. Adding 0 m^{2}s^{-1} to Diff_coeffs array.')
    settings.Diff_coeffs = [0,settings.Diff_coeffs];
end

if ~any(isnan(settings.Noise))
    disp('We require to simulate zero noise for lookup table calculation. Adding no noise to Noise array.')
    settings.Noise = [NaN,settings.Noise];
end

% Check previous lookup table exists
if settings.Use_Previous_Lookup == 1
    if ~exist([settings.filepath,filesep,settings.lookup_filename])
        settings.Use_Previous_Lookup = 0;
        disp('Use_Previous_Lookup is true, but cannot find a lookup table. Switching to false.')
    end
end

if settings.HR_TR == 1 && settings.Use_Previous_Lookup == 0
    disp('Use_Previous_Lookup is false, but HR TR is turned on. Switching to true.')
    settings.Use_Previous_Lookup = 1;
end

if isfield(settings,'Additional_Loop_Counter') && settings.Use_Previous_Lookup == 0
    %disp('Additional loop counter is on, but Use_Previous_Lookup is false. Switching Use_Previous_Lookup to true.')
    %settings.Use_Previous_Lookup = 1;
    warning('Additional loop counter is on, but Use_Previous_Lookup is false.')
end

if ~exist(settings.filepath,'dir')
    mkdir(settings.filepath);
end

% Optional - Overwrite sequence specific settings for additional looping of parameters
% e.g. useful for testing different sequence TRs etc.
if isfield(settings,'LoopValues') && isfield(settings,'Additional_Loop_Counter')
    settings.(settings.LoopFieldName) = settings.LoopValues(settings.Additional_Loop_Counter,:);
end

% Calculate slice ordering
settings = Calc_Slice_Order(settings);

% Check if simulations already ran
[already_ran,settings.savefilename] = check_if_already_ran(settings);

if already_ran % Don't re-simulate if input parameters unchanged
    disp('Input parameters unchanged - results not re-simulated.')
    load(fullfile(settings.filepath,settings.savefilename),'results','settings')
else % Simulate sequence
    disp('Starting simulations.')
    [simulation_results,settings] = seq_loop(settings); 
    disp('Simulations successful!')
    
    % Process simulation results
    disp('Starting analysis.')
    [results] = analysis_function_core(simulation_results, settings,results);
    disp('Analysis finished.')
    
    if length(settings.Dynamic_Range) > 100 && settings.UseSyntheticData == 0
        [Dynamic_Range,DR_Values] = Calc_Dynamic_Range(results,settings,plot_settings);
        results.Dynamic_Range = Dynamic_Range;
        results.DR_Values = DR_Values;
    end
    
    if isfield(settings,'Additional_Loop_Counter')
        settings.savefilename = [settings.filename,'_Loop',num2str(settings.Additional_Loop_Counter)];
    end
    save(fullfile(settings.filepath,[settings.savefilename,'.mat']),'results','settings')
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


