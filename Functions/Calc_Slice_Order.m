function settings = Calc_Slice_Order(settings)
% Calculate the slice acquisition order

if strcmp(settings.Slice_Order_Type,'OddThenEven')
    settings.Slice_Order = [1:2:settings.Scan_Size(2),2:2:settings.Scan_Size(2)];
elseif strcmp(settings.Slice_Order_Type,'EvenThenOdd')
    settings.Slice_Order = [2:2:settings.Scan_Size(2),1:2:settings.Scan_Size(2)];
elseif strcmp(settings.Slice_Order_Type,'Linear')
    settings.Slice_Order = 1:settings.Scan_Size(2);  
elseif strcmp(settings.Slice_Order_Type,'CentreOut')
    A = (settings.Centre_Slice+1):1:settings.Scan_Size(2);
    B = (settings.Centre_Slice-1):-1:1;
    N = min(numel(A),numel(B));
    settings.Slice_Order = [settings.Centre_Slice,reshape([B(1:N);A(1:N)],1,[]),A(N+1:end),B(N+1:end)];
elseif strcmp(settings.Slice_Order_Type,'CentreIn')
    
else
    settings.Slice_Order = 'NA';
end


end

