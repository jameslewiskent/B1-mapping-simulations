function [results,settings] = run_sequence_simulations(settings,results,plot_settings)
% Function handles iterating through scheme(s) to simulate and hands over to analysing
% data in analysis_function_core.
%
% James Kent. 2023. Using B. Hargreaves EPG Functions.
% TODO: put warnings/verbose info into seperate function instead of being
% here (or in sequence settings)

settings.filepath = fullfile('Data',lower(settings.Scheme));
settings = Create_Filename(settings); % Generate filename

% Load sequence specific settings
settings = load_sequence_settings(settings);

% Optional - Overwrite sequence specific settings for additional looping of parameters
% e.g. useful for testing different sequence TRs etc.
if isfield(settings,'LoopValues') && isfield(settings,'Additional_Loop_Counter')
    settings.(settings.LoopFieldName) = settings.LoopValues(settings.Additional_Loop_Counter,:);
end
if isfield(settings,'LoopValues2') && isfield(settings,'Additional_Loop_Counter2')
    settings.(settings.LoopFieldName2) = settings.LoopValues2(settings.Additional_Loop_Counter2,:);
end

settings.Scan_Size(1) = round(settings.PE1_Resolution*settings.PE1_Partial_Fourier*settings.Matrix_Size(1));
settings.Scan_Size(2) = round(settings.PE2_Resolution*settings.PE2_Partial_Fourier*settings.Matrix_Size(2));
settings = Calc_Slice_Positions(settings); % Calculate Slice positions

if strcmpi(settings.MSor3D,'3D')
    settings.Slice_Shifts = zeros(1,settings.Scan_Size(2));
    settings.PP_Shifts = zeros(1,settings.Scan_Size(2));
elseif strcmpi(settings.MSor3D,'2D')
    settings.PP_Shifts = settings.Slice_Shifts;
end

% Calculate number of TRs
settings.N_TRs = (settings.Dummy_Scans + settings.Segment_Factor*settings.Scan_Size(2))*settings.Modes; % Sandwich & SA2RAGE
if strcmpi(settings.Scheme,'SatTFL')
    settings.N_TRs = 2*settings.N_TRs; % SatTFL
end

settings.prep_spoils = 1.* ones(1,settings.N_TRs); % Number of unit gradients to move through after pre-pulse (spoiling- need enough for optional dummy scans too)
settings.train_spoils = 1.* ones(1,ceil(settings.Scan_Size(1)/settings.Segment_Factor));

% Calculate segment sizes
settings = Calc_Segment_Sizes(settings);

if strcmpi(settings.Scheme,'DREAM')
    settings = Calc_Slice_Delay_Time(settings); % Need to calculate time between neighbouring slices
else
    settings = Check_Min_TR(settings); % Check minimum TR time
end

if ~(strcmpi(settings.UseSyntheticData,'Duke') || strcmpi(settings.UseSyntheticData,'Phantom'))
    settings.Mag_Track_Flags = zeros(1,size(settings.Mag_Track_FAValues,2));
elseif (strcmpi(settings.UseSyntheticData,'Duke') || strcmpi(settings.UseSyntheticData,'Phantom'))
    settings.Mag_Track_Flags = 0;
end
settings.lookup_filename = [settings.Scheme,'_lookup_table.mat'];

if (strcmpi(settings.UseSyntheticData,'Duke') || strcmpi(settings.UseSyntheticData,'Phantom')) && settings.Modes < 8  && strcmp(settings.Enc_Scheme,'Indiv')
    disp('Multitransmit mapping is on, and encoding scheme is individual but the number of modes is less than the number of channels (8). Setting modes to 8.')
    settings.Modes = 8;
elseif (strcmpi(settings.UseSyntheticData,'Duke') || strcmpi(settings.UseSyntheticData,'Phantom')) && settings.Modes > 8  && strcmp(settings.Enc_Scheme,'Indiv')
    disp('Multitransmit mapping is on, and encoding scheme is individual but the number of modes is more than the number of channels (8). Setting modes to 8.')
    settings.Modes = 8;
end

if (strcmpi(settings.UseSyntheticData,'Duke') || strcmpi(settings.UseSyntheticData,'Phantom'))
    settings.Mag_Track_FAValues = 1;
end

if (strcmpi(settings.UseSyntheticData,'Duke') || strcmpi(settings.UseSyntheticData,'Phantom')) && settings.Modes == 1 && ~strcmp(settings.Enc_Scheme,'CP')
    disp('For Modes = 1 encoding scheme is fixed as CP. Setting settings.Enc_Scheme = "CP"')
    settings.Enc_Scheme = 'CP';
end

if (strcmpi(settings.UseSyntheticData,'Duke') || strcmpi(settings.UseSyntheticData,'Phantom')) && settings.Repeats > 1
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

if isfield(settings,'Additional_Loop_Counter') && settings.Use_Previous_Lookup == 0 && (strcmpi(settings.Scheme,'Sandwich') || strcmpi(settings.Scheme,'SA2RAGE'))
    %disp('Additional loop counter is on, but Use_Previous_Lookup is false. Switching Use_Previous_Lookup to true.')
    %settings.Use_Previous_Lookup = 1;
    warning('Additional loop counter is on, but Use_Previous_Lookup is false.')
end

if ~exist(settings.filepath,'dir')
    mkdir(settings.filepath);
end

if settings.Global_T1 == 1 && (strcmpi(settings.UseSyntheticData,'Duke') || strcmpi(settings.UseSyntheticData,'Phantom'))
    % Overwrite synthetic T1s with a single global value
    disp('Fixed global T1 value.');
elseif (strcmpi(settings.UseSyntheticData,'Duke') || strcmpi(settings.UseSyntheticData,'Phantom'))
    settings.T1s = 0;
    disp('Using synthetic T1s, not those specified in settings.T1s.');
end

if strcmpi(settings.Scheme,'GRE') && ~strcmpi(settings.LoopFieldName2,'Coil_Cycle_Order')
    if strcmpi(settings.Coil_Cycle,'SS') || strcmpi(settings.Coil_Cycle,'SW')
        settings.Coil_Cycle_Order = 1:8;
    elseif strcmpi(settings.Coil_Cycle,'CC')
        settings.Coil_Cycle_Order = [1,4,7,2,5,8,3,6];
    end
end

if isfield(settings,'Coil_Cycle') && isfield(settings,'Coil_Cycle_Order')
    if strcmpi(settings.Coil_Cycle,'CC')
        disp(['Coil cycling is turned on. Coil-order: ',num2str(settings.Coil_Cycle_Order)]);
    else
        disp(['Coil cycling is turned off. Coil-order: ',num2str(settings.Coil_Cycle_Order)]);
    end
end

if strcmpi(settings.Scheme,'GRE')
    % Must get RF pulse after additional loop counters, incase e.g. FA etc. is changed
    % TODO: Need to update this for other sequences too
    settings.RF_Pulse = Get_RF_Pulse(settings.nom_FA,settings.RF_Type,settings.RF_Time,settings.RF_TBP,settings.Ref_Voltage);
end

if (strcmpi(settings.UseSyntheticData,'Duke') || strcmpi(settings.UseSyntheticData,'Phantom'))
    % Now pulse voltages are known, we can calculate the associated transmit field
    
    % These are the ground truth mode transmit maps
    [settings.Dynamic_Range,settings.Tx_FA_map,settings.Enc_Mat] = Calc_Tx(max(settings.RF_Pulse),settings);
    disp(['Multi-transmit mode mapping is active. Simulating B1 Mapping of ', num2str(size(settings.Tx_FA_map,3)),' Transmit Mode Configuration.']);
    
    % These are the ground truth single-channel maps
    if settings.Modes > 1 && ~strcmpi(settings.Enc_Scheme,'Indiv')
        [~,settings.GT_Tx_FA_map,~] = Calc_Tx(max(settings.RF_Pulse),settings,Calc_Enc_Mat('Indiv',settings.Modes));
    else
        settings.GT_Tx_FA_map = settings.Tx_FA_map;
    end
end

% Calculate slice ordering
settings = Calc_Slice_Order(settings);

% Check if simulations already ran
[already_ran,settings.filename] = check_if_already_ran(settings);

if already_ran % Don't re-simulate if input parameters unchanged
    disp('Input parameters unchanged - results not re-simulated.')
    load(fullfile(settings.filepath,settings.filename),'results')
    
else % Simulate sequence
    [simulation_results,settings] = seq_loop(settings);
    
    % Process simulation results
    [results] = analysis_function_core(simulation_results, settings,results);
    
    if length(settings.Dynamic_Range) > 100 && ~(strcmpi(settings.UseSyntheticData,'Duke') || strcmpi(settings.UseSyntheticData,'Phantom'))
        [Dynamic_Range,DR_Values] = Calc_Dynamic_Range(results,settings,plot_settings);
        results.Dynamic_Range = Dynamic_Range;
        results.DR_Values = DR_Values;
    end
    
    temp_settings = settings; % Create temporary store of settings before removing some fields
    try
        % Remove loop field name and values from saved settings because this is irrelevant and
        % prevents additional loop information being overwritten in case
        % previous results/settings are loaded in
        settings = rmfield(settings,'LoopFieldName');
        settings = rmfield(settings,'LoopValues');
        settings = rmfield(settings,'Additional_Loop_Counter');
    end
    try
        % Also remove loop field name and values because this is irrelevant
        settings = rmfield(settings,'LoopFieldName2');
        settings = rmfield(settings,'LoopValues2');
        settings = rmfield(settings,'Additional_Loop_Counter2');
    end
    
    % Save results
    save(fullfile(settings.filepath,[settings.filename,'.mat']),'results','settings')
    settings = temp_settings; % Restore removed fields
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


