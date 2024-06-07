function settings = Calc_Slice_Delay_Time(settings)
% Need to calculate time between neighbouring slices

% Calculate slice ordering
settings = Calc_Slice_Order(settings);
Num_Slices_Between_Adjacent = find(settings.Slice_Order == 2) - find(settings.Slice_Order == 1) -1; % Find how many slices are between neighbouring slices
if isempty(Num_Slices_Between_Adjacent)
    Num_Slices_Between_Adjacent = 0;
end
Slice_Duration = (settings.Segment_Sizes(1)*settings.IT_TR + settings.Td); % Minimum duration of slice acquisition

% Delay time split across all slices
settings.Compulsory_Delay_Time = max(0,(settings.Neighbouring_Slice_Delay_Time - (Num_Slices_Between_Adjacent*Slice_Duration))/(Num_Slices_Between_Adjacent+1));

end

