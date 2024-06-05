function str = num2strpn(num,val)
% Displays postive or negative numbers as string with associated sign
if num > 0
    str = ['+ ',num2str(num,val)];
else
    str = ['- ',num2str(abs(num),val)];
end

end