function [x_ss] = Calc_SS(y)
% Determine if and when a steady state is reached
y = diff(y);
indices = find(y<0);

y = y(indices);
y = y(y<-0.1);

y = abs([diff(y),0]./y);

SS_index = find(y<0.01,1);
x_ss = indices(SS_index);
end

