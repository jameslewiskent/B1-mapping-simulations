function settings = Check_Mag_Track_Flag(Dynamic_Range_n,T1_n,B0_n,Flow_n,Diff_n,settings)
% Check if magnetisation tracking should be enabled

% Range of FA values
if ~strcmpi(settings.Scheme,'AFI')
    FARange = settings.Dynamic_Range(1,:)*settings.nomPP_FA*(180/pi);
elseif strcmpi(settings.Scheme,'AFI')
    FARange = settings.Dynamic_Range(1,:)*settings.nomFA*(180/pi);
end

% Find index of DR that is closest to the requested FA values for
% magnetisation tracking
Mag_Track_DRValues = zeros(1,size(settings.Mag_Track_FAValues,2));
for FA_n = 1:size(settings.Mag_Track_FAValues,2)
    [min_val,Mag_Track_DRValues(FA_n)] = min(abs(settings.Mag_Track_FAValues(FA_n) - FARange));
    if min_val > 10
        Mag_Track_DRValues(FA_n) = nan;
    end
end



if any(Mag_Track_DRValues == Dynamic_Range_n) && any(settings.Mag_Track_T1Values == settings.T1s(T1_n)) && settings.B0_Range_Hz(B0_n) == 0 && settings.Velocities(Flow_n) == 0 && settings.Diff_coeffs(Diff_n) == 0
    settings.Mag_Track_Flags = Mag_Track_DRValues == Dynamic_Range_n;
    %disp('Tracking Magnetisation')
elseif settings.UseSyntheticData == 1
    % Convert given index to value in long array
    Mag_Track_Long_SynInd = sub2ind([139,178],settings.Mag_Track_SynInd(1),settings.Mag_Track_SynInd(2));
    if Dynamic_Range_n == Mag_Track_Long_SynInd
        settings.Mag_Track_Flags = 1;
    end
else
    settings.Mag_Track_Flags = zeros(1,size(settings.Mag_Track_FAValues,2));
end

end

