function Generate_Lookup_Table(Image_Train_Ratio,settings)
% Generate a lookup table from simulation, then apply that lookup table to
% the rest of the data

Image_Train_Ratio = real(Image_Train_Ratio);

% if isfile([settings.filepath,filesep,settings.lookup_filename])
%     load([settings.filepath,filesep,settings.lookup_filename],'x_query','fx_interp');
% end
% if isfile([settings.filepath,filesep,settings.lookup_filename]) && size(x_query,2) == settings.Lookup_Size && 
%     % If lookup table exists and the size requested is what is expected,
%     % then we can use the existing table
%     disp('Using a previously saved lookup table.')
% else
%     disp('No lookup table was found. Generating new lookup table.')
    disp('Generating lookup table.')
    
    % Choosen T1, On resonance, No diffusion, flow or noise and an average over repeats. 
    try
    fx = squeeze(mean(Image_Train_Ratio(:,:,:,:,settings.B0_Range_Hz == 0,settings.T1s == settings.Lookup_T1,settings.Velocities == 0,settings.Diff_co == 0,isnan(settings.Noise),:),10));
    
    fx(isnan(fx)) = 1; % Replace nans with 1 (Fixes lookup for FA = 0 divide by 0 values)
    catch
        error('Could not create a lookup table with B0 = 0 Hz, Diff_co = 0, Flow = 0 and Noise = 0. Are we simulating these values?')
    end
    [~,minind] = min(fx); fx(minind:end) = fx(minind); % Ensure monotonic
    
    % Interpolate function
    x_query = linspace(settings.Dynamic_Range(1)*settings.nomPP_FA,settings.Dynamic_Range(end)*settings.nomPP_FA,settings.Lookup_Size); % FA in degrees
    fx_interp = interp1(settings.Dynamic_Range*settings.nomPP_FA, fx, x_query, 'spline');
    
    %figure(); plot(x_query,fx_interp)
    disp(['Lookup table saved: ',[settings.filepath,filesep,settings.lookup_filename]])
    save([settings.filepath,filesep,settings.lookup_filename],'x_query','fx_interp');
    
    generateLookupCPPandHeaderFiles(settings.PP_RF_Type,x_query,fx_interp); % Save lookup table as .cpp and .h files for sequence B1 map reconstruction
end