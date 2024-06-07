function settings = Check_Min_TR(settings)
% Need to calculate scan duration, to ensure a minimum duration between subsequent slices
Reference_Duration = settings.Scan_Size(2)*(settings.Segment_Sizes(1)*settings.IT_TR + settings.PP_RF_Time); % Minimum duration of all reference images
if Reference_Duration > settings.TR
    settings.Compulsory_Delay_Time = settings.Min_Delay_Time;
else
    settings.Compulsory_Delay_Time = (settings.TR - Reference_Duration)/(settings.Scan_Size(2)*settings.Segment_Factor);
end

end

