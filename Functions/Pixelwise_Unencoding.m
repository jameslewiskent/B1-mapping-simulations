function [outputArg1,outputArg2] = Pixelwise_Unencoding(Enc_Mat,W_Mat,Modes)
% Do pixel-by-pixel un-encoding

p = 0.25; % Empirically determined constants for weighting matrix
q = 20; % Empirically determined constants for weighting matrix
Channels = 8; % Channels

for indi = 1:size(slice,1)
    for indj = 1:size(slice,2)
        
        if slice_T1(indi,indj) ~= 0 % only run for simulated data
            for B0_n = 1:length(TRs)
                for T1_n = 1 % sim only ran for one T1 value (from synthetic data)
                    for Noise_n = 1:length(Noise_pcs)
                        for channel = 1:Channels
                            B1_channel_sum = 0;
                            IT = zeros(1,Channels);
                            
                            % Calculate W matrix from pixel value
                            for m1 = 1:Modes
                                for m2 = 1:Modes % due to brace indexing, have to go through each value of cell and pull out values across the modes
                                    IT(m2) = Kappa_error_results{3,m2}(indi,indj,B0_n,T1_n,Noise_n);
                                end
                                Jm = abs(IT(m1)) / max(abs(IT));
                                W_Mat(m1,m1) = 1/(1 + exp(-q*(Jm-p)));
                            end
                            % Now calculate un-encoding matrix for specific pixel
                            Un_Enc_Mat = ((pinv(Enc_Mat'*W_Mat*Enc_Mat))*Enc_Mat'*W_Mat);
                            
                            
                            for mode = 1:Modes
                                %On a voxel by voxel basis, multply
                                %unencoding matrix by the resulting b1
                                %maps (which gets its phase from the
                                %reference image)
                                B1_channel_sum =  B1_channel_sum + Un_Enc_Mat(channel,mode).*Kappa_error_results{1,mode}(indi,indj,TR_n,T1_n,Noise_n).*exp(1i*angle(conj(Kappa_error_results{3,mode}(indi,indj,TR_n,T1_n,Noise_n))));
                            end
                            
                            Kappa_error_results{5,channel}(indi,indj,B0_n,T1_n,Noise_n) = B1_channel_sum ; % Store channel maps in with mode maps
                            
                        end
                    end
                end
            end
            
        else
            Kappa_error_results{5,channel}(indi,indj,:,:,:) = 0 ; % Store channel maps in with mode maps
        end
        
        
    end
end
end

