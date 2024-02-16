function Measured_FA = Apply_Lookup_Table(Image_Train_Ratio,settings)
% Apply that lookup table to the data
Image_Train_Ratio = real(Image_Train_Ratio);

load([settings.filepath,filesep,settings.lookup_filename],'x_query','fx_interp');

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
