function [Max_Val_IT,FWHM] = Process_IT(settings,IT,ITn)
% Process an individual image train

% Re-order Image Trains
% Collected chronologically in the image train.
if strcmpi(settings.Scheme,'DREAM')  % Don't re-order DREAM k-space
    Reordered_IT = IT;
else
    Reordered_IT = zeros(size(IT));

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
        Reordered_IT(Reorder_PE1(ITn,IT_n),:,:,:,:,:,:,:,:,:) = IT(IT_n,:,:,:,:,:,:,:,:,:);
    end
end
clearvars IT Reorder

if settings.UseSyntheticData == 0
    Reordered_IT = Reordered_IT(:,settings.Centre_Slice,:,:,:,:,:,:,:,:); % Select Centre slice
end

% Add zero-mean complex Gaussian noise independently to each reordered IT
%Noise_sd = max(abs(Reordered_IT1),[],'all')./db2mag(settings.Noise);
Noise_sd = 1./db2mag(settings.Noise);
Noise_sd(isnan(settings.Noise)) = 0; % Set any NaN's to 0, no noise added

IT_Noisy = Reordered_IT + bsxfun(@times,randn([size(Reordered_IT,1:8),size(settings.Noise,2),settings.Repeats]),reshape(Noise_sd,[1 1 1 1 1 1 1 1 size(settings.Noise,2) 1])) + 1i.*bsxfun(@times,randn([size(Reordered_IT,1:8),size(settings.Noise,2),settings.Repeats]),reshape(Noise_sd,[1 1 1 1 1 1 1 1 size(settings.Noise,2) 1]));
clearvars Reordered_IT

% Zero-fill PE1, accounts for partial Fourier by unevenly filling
if settings.PE1_Partial_Fourier ~= 1 || settings.PE1_Resolution ~= 1
    Centre_Line = ceil(settings.PE1_Resolution*(settings.Matrix_Size(1)/2 - (settings.Matrix_Size(1)*(1-settings.PE1_Partial_Fourier)))); % (zero-indexed)
    Pad_Size = settings.Matrix_Size(1)-size(IT_Noisy,1);
    PadTop_Size = (settings.Matrix_Size(1)/2)-Centre_Line;
    PadBot_Size = Pad_Size-PadTop_Size;
    PadTop = zeros([PadTop_Size,size(IT_Noisy,2:ndims(IT_Noisy))]);
    PadBot = zeros([PadBot_Size,size(IT_Noisy,2:ndims(IT_Noisy))]);
    IT_Noisy = cat(1,cat(1,PadTop, IT_Noisy),PadBot);
end
clearvars PadTop PadBot

% Zero-fill PE2, accounts for partial Fourier by unevenly filling
if (settings.PE2_Partial_Fourier ~= 1 || settings.PE2_Resolution ~= 1) && settings.Matrix_Size(2) > 1
    Centre_Line = ceil(settings.PE2_Resolution*(settings.Matrix_Size(2)/2 - (settings.Matrix_Size(2)*(1-settings.PE2_Partial_Fourier)))); % (zero-indexed)
    Pad_Size = settings.Matrix_Size(2)-size(IT_Noisy,1);
    PadTop_Size = (settings.Matrix_Size(2)/2)-Centre_Line;
    PadBot_Size = Pad_Size-PadTop_Size;
    PadTop = zeros([size(IT_Noisy,1),PadTop_Size,size(IT_Noisy,3:ndims(IT_Noisy))]);
    PadBot = zeros([size(IT_Noisy,1),PadBot_Size,size(IT_Noisy,3:ndims(IT_Noisy))]);
    IT_Noisy = cat(2,cat(1,PadTop, IT_Noisy),PadBot);
end
clearvars PadTop PadBot

% Fourier transform image train
if strcmp(settings.MSor3D,'3D')
    fftdims = [1 2];
elseif strcmp(settings.MSor3D,'2D')
    fftdims = [1];
end
for f = fftdims
    IT_Noisy =  ifftshift(ifft(fftshift(IT_Noisy,f),[],f),f);
end
FT_IT = IT_Noisy;
clearvars IT_Noisy

if settings.UseSyntheticData == 0 && size(FT_IT,4) ~= 1 && ~isfield(settings,'LoopValues')
    % Evaluate abs FWHM for all readout
    centre_index = floor(size(FT_IT,2)/2 +1);
    figure('Name',['FFT of IT for ',settings.Scheme,' sequence'],'color','w')
    tiles = tiledlayout('flow','tilespacing','compact','padding','none');
    title(tiles,settings.title_string)
    nexttile; mesh(abs(squeeze(FT_IT(:,centre_index,1,:,1,1,1,1,1,1,1)))) % try mesh or imagesc
    caxis([0 7])
    title('Image Train 1')
    ylabel('Pre-pulse FA');
    xlabel('Index');
end

FWHM = zeros([2,2,size(FT_IT,4:ndims(FT_IT))]);
if settings.Calc_FWHM == 1 % Measure FWHM (can take a while)
    for Dynamic_Range_n = 1:size(FT_IT,4)
        for B0_Range_n = 1:size(FT_IT,5)
            for T1_n = 1:size(FT_IT,6)
                for Flow_n = 1:size(FT_IT,7)
                    for Diff_n = 1:size(FT_IT,8)
                        for Noise_n = 1:size(FT_IT,9)
                            for Repeat_n = 1:size(FT_IT,10)
                                FWHM(1,1,Dynamic_Range_n,B0_Range_n,T1_n,Flow_n,Diff_n,Noise_n,Repeat_n) = findFWHM(FT_IT(:,floor(size(FT_IT,2)/2 +1),:,Dynamic_Range_n,B0_Range_n,T1_n,Flow_n,Diff_n,Noise_n,Repeat_n));
                                if size(FT_IT,2) > 1
                                    FWHM(1,2,Dynamic_Range_n,B0_Range_n,T1_n,Flow_n,Diff_n,Noise_n,Repeat_n) = findFWHM(FT_IT(floor(size(FT_IT,1)/2 +1),:,:,Dynamic_Range_n,B0_Range_n,T1_n,Flow_n,Diff_n,Noise_n,Repeat_n)');
                                end
                            end
                        end
                    end
                end
            end
        end
        if settings.verbose == 1
            disp(['Calculated FWHM for DR: ',num2str(Dynamic_Range_n),' of ',num2str(size(FT_IT,4)), ' completed.'])
        end
    end
end

if settings.Sum_PSF == 0 
    % Take centre of PSF
    Max_Val_IT = max(FT_IT(floor(size(FT_IT,1)/2 +1),floor(size(FT_IT,2)/2 +1),:,:,:,:,:,:,:,:),[],1); % Find index of maximum value
elseif settings.Sum_PSF == 1
    % Sum PSF
    Max_Val_IT = sum(FT_IT,1);
end
clearvars FT_IT index
end

