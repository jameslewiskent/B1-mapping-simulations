function FWHM = findFWHM(fx)
% Finds the full width half maximum of function

fx = squeeze(abs(fx'));
dx = 1e-2;
x = 1:size(fx,2);
x_query = 1:dx:size(x,2);

[m, ~] = max(fx);		%	Find maximum value and index
fx_interp = interp1(x, fx, x_query, 'spline'); % interpolate function

ind = find(fx_interp>=m/2);	%	Find indices just above max(fx)/2
nl = min(ind);			%	Leftmost index
nr = max(ind);			%	Rightmost index

%	Get FWHM
FWHM = x_query(nr)-x_query(nl);


% figure(); plot(x,fx,'ro'); hold on
% plot(x_query,fx_interp,'b');
% plot(x_query(nl),fx_interp(nl),'gx');
% plot(x_query(nr),fx_interp(nr),'gx');

end