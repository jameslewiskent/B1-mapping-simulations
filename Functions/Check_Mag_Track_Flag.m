function settings = Check_Mag_Track_Flag(Dynamic_Range_n,T1_n,B0_n,settings)
% Check if magnetisation tracking should be enabled

% Range of FA values
if ~strcmp(settings.Scheme,'AFI')
FARange = settings.Dynamic_Range*settings.nomPP_FA*(180/pi);
elseif strcmp(settings.Scheme,'AFI')
FARange = settings.Dynamic_Range*settings.nomFA*(180/pi);
end

% Find index of DR that is closest to the requested FA values for
% magnetisation tracking
Mag_Track_DRValues = zeros(1,size(settings.Mag_Track_FAValues,2));
for FA_n = 1:size(settings.Mag_Track_FAValues,2)
    [~,Mag_Track_DRValues(FA_n)] = min(abs(settings.Mag_Track_FAValues(FA_n) - FARange));
end

if any(Mag_Track_DRValues == Dynamic_Range_n) && any(settings.Mag_Track_T1Values == settings.T1s(T1_n)) && settings.B0_Range_Hz(B0_n) == 0
settings.Mag_Track_Flags = Mag_Track_DRValues == Dynamic_Range_n;
disp('Tracking Magnetisation')
else
settings.Mag_Track_Flags = zeros(1,size(settings.Mag_Track_FAValues,2));
end

end

