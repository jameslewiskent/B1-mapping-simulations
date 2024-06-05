function settings = Calc_Slice_Positions(settings)
% Calculate the slice positions (cm)

settings.Centre_Slice = floor((settings.Scan_Size(2)/2) +1);
settings.Slice_Positions = zeros(settings.Scan_Size(2),3);
settings.Slice_Positions(:,1) = ((1:settings.Scan_Size(2))-settings.Centre_Slice)*settings.Slice_Thickness;
end