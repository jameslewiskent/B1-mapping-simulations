function settings = Calc_Slice_Order(settings)
% Calculate the slice acquisition order

if strcmp(settings.Slice_Order_Type,'OddThenEven')
    settings.Slice_Order = [1:2:settings.Scan_Size(2),2:2:settings.Scan_Size(2)];
elseif strcmp(settings.Slice_Order_Type,'EvenThenOdd')
    settings.Slice_Order = [2:2:settings.Scan_Size(2),1:2:settings.Scan_Size(2)];
elseif strcmp(settings.Slice_Order_Type,'Linear')
    settings.Slice_Order = 1:settings.Scan_Size(2);  
else
    settings.Slice_Order = 'NA';
end


end

