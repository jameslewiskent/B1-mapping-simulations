function Measured_FA = Lookup_Table(Image_Train_Ratio,settings)

Image_Train_Ratio = real(Image_Train_Ratio);

% Choose T1, On resonance and average over repeats. Choose no diffusion or flow
fx = squeeze(mean(Image_Train_Ratio(:,:,:,:,settings.B0_Range_Hz == 0,settings.T1s == settings.Lookup_T1,settings.Velocities == 0,settings.Diff_co == 0,1,:),10));

[~,minind] = min(fx); fx(minind:end) = fx(minind); % Ensure monotonic

% Interpolate function
x_query = linspace(settings.Dynamic_Range(1),settings.Dynamic_Range(end),1e3);
fx_interp = interp1(settings.Dynamic_Range, real(fx), x_query, 'spline');

Measured_FA = zeros(size(Image_Train_Ratio));
for Mode_n = 1:size(settings.Tx_FA_map,3)
    for Dynamic_Range_n = 1:size(Image_Train_Ratio,4)
        for B0_n = 1:size(Image_Train_Ratio,5)
            for T1_n = 1:size(Image_Train_Ratio,6)
                for Flow_n = 1:size(Image_Train_Ratio,7)
                    for Diff_n = 1:size(Image_Train_Ratio,8)
                        for Noise_n = 1:size(Image_Train_Ratio,9)
                            for Repeat_n = 1:size(Image_Train_Ratio,10)
                                [~,min_ind] = min(abs(Image_Train_Ratio(1,1,1,Dynamic_Range_n,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Repeat_n) - fx_interp));
                                Measured_FA(1,1,Mode_n,Dynamic_Range_n,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Repeat_n) = x_query(min_ind)*settings.nomPP_FA; % Measured_FA in radians
                            end
                        end
                    end
                end
            end
        end
    end
end

end
