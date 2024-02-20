function [Unencoded_Image_Maps] = Pixelwise_Unencoding(settings,Image_Maps)
% JK Do pixel-by-pixel un-encoding following DREAM paper
% Image_Maps input should already be complex
disp('Unencoding maps.')

p = 0.25; % Empirically determined constants for weighting matrix
q = 20; % Empirically determined constants for weighting matrix
Channels = 8; % Channels (Fixed from synthetic data)

Enc_Mat = settings.Enc_Mat;
Unencoded_Image_Maps = zeros([size(Image_Maps,1:2),Channels,size(Image_Maps,4:ndims(Image_Maps))]);
for indi = 1:size(Image_Maps,1)
    for indj = 1:size(Image_Maps,2)
        if settings.Synthetic_Mask(indi,indj)
            for B0_n = 1:length(settings.B0_Range_Hz)
                for T1_n = 1:length(settings.T1s)
                    for Flow_n = 1:length(settings.Velocities)
                        for Diff_n = 1:length(settings.Diff_coeffs)
                            for Noise_n = 1:length(settings.Noise)
                                for Repeat_n = 1:length(settings.Repeats)
                                    for Channel_n = 1:Channels
                                        if settings.Modes > Channels
                                            % Calculate W matrix from pixel value
                                            for Mode_n = 1:settings.Modes
                                                IT = squeeze(Image_Maps(indi,indj,:,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Repeat_n));
                                                Jm = abs(IT(Mode_n)) / max(abs(IT));
                                                W_Mat(Mode_n,Mode_n) = 1/(1 + exp(-q*(Jm-p)));
                                            end
                                        else
                                            W_Mat = eye(settings.Modes,settings.Modes);
                                        end
                                        % Now calculate un-encoding matrix for specific pixel
                                        Un_Enc_Mat = ((pinv(Enc_Mat'*W_Mat*Enc_Mat))*Enc_Mat'*W_Mat);
                                        
                                        B1_channel_sum = 0;
                                        for Mode_n = 1:settings.Modes
                                            % On a voxel by voxel basis, multiply unencoding matrix by the resulting b1 maps (which has already got its phase from the reference image)
                                            B1_channel_sum =  B1_channel_sum + Un_Enc_Mat(Channel_n,Mode_n)'.*Image_Maps(indi,indj,Mode_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Repeat_n);
                                        end
                                        Unencoded_Image_Maps(indi,indj,Channel_n,1,B0_n,T1_n,Flow_n,Diff_n,Noise_n,Repeat_n) = B1_channel_sum ; % Store channel maps in with mode maps
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
disp('Finished un-encoding maps.')

end

