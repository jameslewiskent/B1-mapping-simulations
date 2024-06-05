function [Body_Mask, Heart_Mask, Synthetic_T1s] = Create_Synthetic_Masks(SyntheticDuke,Slice_Number,settings)
% Create masks and T1 maps of synthetic body data.
% Assign T1 values to synthetic data

slice = SyntheticDuke.sigma(:,:,Slice_Number); % For now only consider one slice of conductivity data to mask

% Define Heart Mask
mask_heart = zeros(size(slice)); % pre-allocate heart mask
mask_blood = zeros(size(slice)); % pre-allocate T2 value array
for indi = 1:size(slice,1)
    for indj = 1:size(slice,2)
        if slice(indi,indj) > 0.85 && slice(indi,indj) < 0.95 % Assume Myocardium
            mask_heart(indi,indj) = 1;
        elseif slice(indi,indj) > 1.3 && slice(indi,indj) < 1.4
            mask_blood(indi,indj) = 1;
        end
    end
end
clear i j
Heart_T1s = mask_heart.*settings.T1_heart + mask_blood.*settings.T1_blood;
Heart_Mask = logical(mask_heart + mask_blood);

% Define whole-body mask regions
mask = zeros(size(slice));
Body_T1s = 2.*double(slice); % Define slice T1 values as 2 * conductivity (very rough)
Body_T1s(round(slice,5,'significant') == round(slice(86,49),5,'significant')) = 0; % Remove Lungs
Body_T1s = Body_T1s.*SyntheticDuke.mask(:,:,Slice_Number); % Remove Tx/Rx

mask(Body_T1s ~= 0) = 1;
Body_Mask = logical(mask);

% Replace T1s less than 0.5 s
Body_T1s(Body_T1s < 0.5 & Body_T1s > 0) = 0.5;

% Replace T1s in heart with specified values
Body_T1s(Heart_Mask) = Heart_T1s(Heart_Mask);
Synthetic_T1s = Body_T1s;

end

