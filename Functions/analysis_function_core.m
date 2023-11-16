function [results] = analysis_function_core(simulation_results, settings,results)
% Function handles the actually analysis i.e. Reordering, zero-filling,
% generating synthetic noise, Fourier transform, calculating alpha.
%
% Only processes one set of image trains at a time (E.g. IT1 and IT2
% for a specific scheme and for multiTx it is ran for however modes there are)

IT1 = simulation_results.IT1;
IT2 = simulation_results.IT2;

% Re-order Image Trains
% Collected chronologically in the image train.
if strcmpi(settings.Scheme,'DREAM')  % Don't re-order DREAM k-space
    Reordered_IT1 = IT1;
    Reordered_IT2 = IT2;
else
    Reordered_IT1 = zeros(size(IT1));
    Reordered_IT2 = zeros(size(IT1));
    
    [Reorder_PE1,Reorder_PE2] = Calc_Reordering(settings);
    
    for Seg_n = 1:size(settings.Segment_Sizes,2)
        % Reorder taking into account any segmentation that has occured.
        Reordered_Seg(1,sum(settings.Segment_Sizes(1:Seg_n - 1))+1 : sum(settings.Segment_Sizes(1:Seg_n))) = Reorder_PE1(1,Seg_n:size(settings.Segment_Sizes,2):end);
        Reordered_Seg(2,sum(settings.Segment_Sizes(1:Seg_n - 1))+1 : sum(settings.Segment_Sizes(1:Seg_n))) = Reorder_PE1(2,Seg_n:size(settings.Segment_Sizes,2):end);
    end
    if size(Reorder_PE1,2) ~= size(Reordered_Seg,2)
        error('ERROR: Segment reordering failed.')
    else
        Reorder_PE1 = Reordered_Seg;
    end
    clearvars Reordered_Seg A B N
    
    for IT_n = 1:settings.PE1_Resolution*settings.PE1_Partial_Fourier*settings.Matrix_Size(1)
        % Perform the reordering of the data
        Reordered_IT1(Reorder_PE1(1,IT_n),:,:,:,:,:,:,:) = IT1(IT_n,:,:,:,:,:,:,:);
        Reordered_IT2(Reorder_PE1(2,IT_n),:,:,:,:,:,:,:) = IT2(IT_n,:,:,:,:,:,:,:);
    end
end
clearvars IT1 IT2 Reorder

if settings.UseSyntheticData == 0
    % Select Centre slice
    Reordered_IT1 = Reordered_IT1(:,settings.Centre_Slice,:,:,:,:,:,:);
    Reordered_IT2 = Reordered_IT2(:,settings.Centre_Slice,:,:,:,:,:,:);
end

% Add zero-mean complex Gaussian noise independently to each reordered IT
%Noise_sd = max(abs(Reordered_IT1),[],'all')./db2mag(settings.Noise);
Noise_sd = 1./db2mag(settings.Noise);
Noise_sd(isnan(settings.Noise)) = 0; % Set any NaN's to 0, no noise added

IT1_Noisy = Reordered_IT1 + bsxfun(@times,randn([size(Reordered_IT1),size(settings.Noise,2),settings.Repeats]),reshape(Noise_sd,[1 1 1 1 1 1 1 1 size(settings.Noise,2) 1])) + 1i.*bsxfun(@times,randn([size(Reordered_IT1),size(settings.Noise,2),settings.Repeats]),reshape(Noise_sd,[1 1 1 1 1 1 1 1 size(settings.Noise,2) 1]));
IT2_Noisy = Reordered_IT2 + bsxfun(@times,randn([size(Reordered_IT2),size(settings.Noise,2),settings.Repeats]),reshape(Noise_sd,[1 1 1 1 1 1 1 1 size(settings.Noise,2) 1])) + 1i.*bsxfun(@times,randn([size(Reordered_IT2),size(settings.Noise,2),settings.Repeats]),reshape(Noise_sd,[1 1 1 1 1 1 1 1 size(settings.Noise,2) 1]));
clearvars Reordered_IT1 Reordered_IT2

% Zero-fill PE1, accounts for partial Fourier by unevenly filling
if settings.PE1_Partial_Fourier ~= 1 || settings.PE1_Resolution ~= 1
    Centre_Line = ceil(settings.PE1_Resolution*(settings.Matrix_Size(1)/2 - (settings.Matrix_Size(1)*(1-settings.PE1_Partial_Fourier)))); % (zero-indexed)
    Pad_Size = settings.Matrix_Size(1)-size(IT1_Noisy,1);
    PadTop_Size = (settings.Matrix_Size(1)/2)-Centre_Line;
    PadBot_Size = Pad_Size-PadTop_Size;
    PadTop = zeros([PadTop_Size,size(IT1_Noisy,2:ndims(IT1_Noisy))]);
    PadBot = zeros([PadBot_Size,size(IT1_Noisy,2:ndims(IT1_Noisy))]);
    IT1_Noisy = cat(1,cat(1,PadTop, IT1_Noisy),PadBot);
    IT2_Noisy = cat(1,cat(1,PadTop, IT2_Noisy),PadBot);
end
clearvars PadTop PadBot

% Zero-fill PE2, accounts for partial Fourier by unevenly filling
if (settings.PE2_Partial_Fourier ~= 1 || settings.PE2_Resolution ~= 1) && settings.Matrix_Size(2) > 1
    Centre_Line = ceil(settings.PE2_Resolution*(settings.Matrix_Size(2)/2 - (settings.Matrix_Size(2)*(1-settings.PE2_Partial_Fourier)))); % (zero-indexed)
    Pad_Size = settings.Matrix_Size(2)-size(IT1_Noisy,1);
    PadTop_Size = (settings.Matrix_Size(2)/2)-Centre_Line;
    PadBot_Size = Pad_Size-PadTop_Size;
    PadTop = zeros([size(IT1_Noisy,1),PadTop_Size,size(IT1_Noisy,3:ndims(IT1_Noisy))]);
    PadBot = zeros([size(IT1_Noisy,1),PadBot_Size,size(IT1_Noisy,3:ndims(IT1_Noisy))]);
    IT1_Noisy = cat(2,cat(1,PadTop, IT1_Noisy),PadBot);
    IT2_Noisy = cat(2,cat(1,PadTop, IT2_Noisy),PadBot);
end
clearvars PadTop PadBot


% Fourier transform image train
for f = [1 2]
    IT1_Noisy =  ifftshift(ifft(fftshift(IT1_Noisy,f),[],f),f);
    IT2_Noisy =  ifftshift(ifft(fftshift(IT2_Noisy,f),[],f),f);
end
FT_IT1 = IT1_Noisy;
FT_IT2 = IT2_Noisy;
clearvars IT1_Noisy IT2_Noisy


if settings.UseSyntheticData == 0 && size(FT_IT1,4) ~= 1
    % Evaluate abs FWHM for all readout
    centre_index = floor(size(FT_IT1,2)/2 +1);
    figure('Name',['FFT of IT for ',settings.Scheme,' sequence'],'color','w')
    tiles = tiledlayout(1,2,'tilespacing','compact','padding','none');
    title(tiles,settings.title_string)
    nexttile; mesh(abs(squeeze(FT_IT1(:,centre_index,1,:,1,1,1,1,1,1,1)))) % try mesh or imagesc
    caxis([0 7])
    title('Image Train 1')
    ylabel('Pre-pulse FA');
    xlabel('Index');
    nexttile; mesh(abs(squeeze(FT_IT2(:,centre_index,1,:,1,1,1,1,1,1,1))))
    caxis([0 7])
    title('Image Train 2')
    ylabel('Pre-pulse FA');
    xlabel('Index');
end

% Measure FWHM
FWHM = zeros([2,size(FT_IT1,4:ndims(FT_IT1))]);
for Dynamic_Range_n = 1:size(FT_IT1,4)
    for B0_Range_n = 1:size(FT_IT1,5)
        for T1_n = 1:size(FT_IT1,6)
            for Flow_n = 1:size(FT_IT1,7)
                for Diff_n = 1:size(FT_IT1,8)
                    for Noise_n = 1:size(FT_IT1,9)
                        for Repeat_n = 1:size(FT_IT1,10)
                            FWHM(1,Dynamic_Range_n,B0_Range_n,T1_n,Flow_n,Diff_n,Noise_n,Repeat_n) = findFWHM(FT_IT1(:,:,:,Dynamic_Range_n,B0_Range_n,T1_n,Flow_n,Diff_n,Noise_n,Repeat_n));
                            FWHM(2,Dynamic_Range_n,B0_Range_n,T1_n,Flow_n,Diff_n,Noise_n,Repeat_n) = findFWHM(FT_IT2(:,:,:,Dynamic_Range_n,B0_Range_n,T1_n,Flow_n,Diff_n,Noise_n,Repeat_n));
                        end
                    end
                end
            end
        end
    end
end

if settings.UseSyntheticData == 0
    % Plot FWHM
    %     figure('Name',['FFT of IT for ',settings.Scheme,' sequence'],'color','w')
    %     tiles = tiledlayout(1,2);
    %     title(tiles,settings.title_string)
    %     nexttile; plot(squeeze(mean(abs(FWHM(1,:,1,1,1,:)),6)),'x') % try mesh or imagesc
    %     caxis([0 7])
    %     title('Image Train 1')
    %     ylabel('Pre-pulse FA');
    %     xlabel('Index');
    %     nexttile; plot(squeeze(mean(abs(FWHM(2,:,1,1,1,:)),6)),'x')
    %     caxis([0 7])
    %     title('Image Train 2')
    %     ylabel('Pre-pulse FA');
    %     xlabel('Index');
end

if settings.Sum_PSF == 0 % Take centre of PSF
    centre_index = floor(settings.Matrix_Size/2 +1);
    Max_Val_IT1 = max(FT_IT1(centre_index(1),:,:,:,:,:,:,:,:,:),[],1); % Find index of maximum value
    Max_Val_IT2 = max(FT_IT2(centre_index(1),:,:,:,:,:,:,:,:,:),[],1); % Store maximum value of IT2
elseif settings.Sum_PSF == 1 % Sum PSF
    Max_Val_IT1 = sum(FT_IT1,1);
    Max_Val_IT2 = sum(FT_IT2,1);
end
clearvars FT_IT1 FT_IT2 index

% Calculate flip angle for various schemes
if strcmpi(settings.Scheme,'DREAM')
    %if strcmp(settings.Echo_Order,'STEFirst')
    Measured_FA = atan(sqrt(abs(2.*(Max_Val_IT1./Max_Val_IT2)))); % Eq. [5] Nehrke, K. et al. MRM. 2012. 68:1517-1526.
    %elseif strcmp(settings.Echo_Order,'FIDFirst')
    %Measured_FA = atan(sqrt(abs(2.*(Max_Val_IT1./Max_Val_IT2)))); % Eq. [5] Nehrke, K. et al. MRM. 2012. 68:1517-1526.
    %end
    
elseif strcmpi(settings.Scheme,'AFI')
    n = settings.TR2/settings.TR1; % TR2/TR1 for AFI Scheme (E.g. n = 100 ms / 20 ms = 5)
    r = Max_Val_IT2./Max_Val_IT1;
    Measured_FA = acos( (r.*n - 1) ./ (n - r) ); % Eq. [6] Yarnykh VL. AFI. Magn Reson Med 2007;57:192?200. https://doi.org/10.1002/mrm.21120.
    
elseif strcmpi(settings.Scheme,'SA2RAGE')
    Measured_FA = Lookup_Table(Max_Val_IT1./Max_Val_IT2,settings);
    
elseif strcmpi(settings.Scheme,'SatTFL')
    if settings.Lookup_T1 ~= 0
        Measured_FA = Lookup_Table(Max_Val_IT2./Max_Val_IT1,settings);
    else
        disp('Calculating flip angle using arcosine instead of lookup table.')
        Measured_FA = acos(Max_Val_IT2./Max_Val_IT1); % Eq. 2 from Chung, S. et al. MRM. 2010. 64(2):439-446.
    end    
elseif  strcmpi(settings.Scheme,'Sandwich')
    if settings.Lookup_T1 ~= 0
        Measured_FA = Lookup_Table(Max_Val_IT2./Max_Val_IT1,settings);
    else
        disp('Calculating flip angle using arcosine instead of lookup table.')
        Measured_FA = acos(Max_Val_IT2./Max_Val_IT1); % Eq. 2 from Chung, S. et al. MRM. 2010. 64(2):439-446.
    end
end

% Pass simulation results
%Max_Val_IT1 Max_Val_IT2 is kept to recon synthetic data mTx maps (using phase of this reference)
results.Measured_FA = (180/pi)*real(Measured_FA);
results.Max_Val_IT1 = Max_Val_IT1;
results.Max_Val_IT2 = Max_Val_IT2;
results.FWHM = FWHM;
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

