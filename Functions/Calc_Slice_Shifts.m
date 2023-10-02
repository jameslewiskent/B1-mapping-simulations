function settings = Calc_Slice_Shifts(settings)
% Calculate the slice shift

if isfield(settings,'Slice_Shift')
disp(['Slice shift value is fixed at ',num2str(settings.Slice_Shift), ' Hz.']);
else
settings.Slice_Shift = settings.Gss*settings.Gamma*settings.Slice_Thickness; % slice shifts [Hz] (G * gamma * dx)
end

settings.Centre_Slice = floor((settings.Scan_Size(2)/2) +1);
settings.Slice_Shifts = ((1:settings.Scan_Size(2))-settings.Centre_Slice)*settings.Slice_Shift;
end

