function settings = Calc_Segment_Sizes(settings)
% Calculate segment sizes

% Calculate number of phase encodes in each segment and check
Train_Size_Tot = settings.Scan_Size(1);
settings.Segment_Sizes = zeros(1,settings.Segment_Factor);
for Seg_n = 1:settings.Segment_Factor
    settings.Segment_Sizes(Seg_n)= ceil(Train_Size_Tot./(settings.Segment_Factor+1-Seg_n));
    Train_Size_Tot =  settings.Scan_Size(1) - sum(settings.Segment_Sizes,'all');
end

% Double-check
if sum(settings.Segment_Sizes,'all') ~= settings.Scan_Size(1)
    error('ERROR: epg_func: Number of phase encodes in segments does not equal size of train requested.')
elseif settings.Segment_Factor ~= 1
    disp(['Number of phase encodes in image train: ',num2str(settings.Scan_Size(1)),', per segment: ',num2str(settings.Segment_Sizes)]);
end

end

