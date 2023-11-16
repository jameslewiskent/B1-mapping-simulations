function Measured_FA = Lookup_Table(Image_Train_Ratio,settings)
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
%end

disp('Applying lookup table.')
Measured_FA = zeros(size(Image_Train_Ratio));
for Mode_n = 1:size(settings.Tx_FA_map,3)
    for Dynamic_Range_n = 1:size(Image_Train_Ratio,4)
        for B0_n = 1:size(Image_Train_Ratio,5)
            for T1_n = 1:size(Image_Train_Ratio,6)
                for Flow_n = 1:size(Image_Train_Ratio,7)
                    for Diff_n = 1:size(Image_Train_Ratio,8)
                        for Noise_n = 1:size(Image_Train_Ratio,9)
                            for Repeat_n = 1:size(Image_Train_Ratio,10)
                                [~,min_ind] = min(abs(Image_Train_Ratio(1,1,Mode_n,Dynamic_Range_n,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Repeat_n) - fx_interp));
                                Measured_FA(1,1,Mode_n,Dynamic_Range_n,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Repeat_n) = x_query(min_ind);
                            end
                        end
                    end
                end
            end
        end
    end
end
disp('Lookup table applied.')
end
