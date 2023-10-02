function [FpFmZ] = epg_trim(FpFmZ,thres)
%
%	Trim higher-order states to N if the sum of absolute
%	values of F+, F- and Z of each order is less than
%	thres for all orders n>N.
%	
%	INPUT:
%		FpFmZ = 3xN vector of F+, F- and Z states.
%		thres = threshold for trimming.
%       OUTPUT:
%               Updated FpFmZ state.
%
%       B.Hargreaves.

f = find(sum(abs(FpFmZ))>=thres);
fn = max(f);

if isempty(f)
    fn = 1; % JK was throwing error when first value below threshold
end

if fn ~= size(FpFmZ,2) && size(FpFmZ,2) >= 10 % Don't trim when there is fewer than 10 states
%disp(['Trimmed from ',num2str(size(FpFmZ,2)),' to ',num2str(fn)])
FpFmZ = FpFmZ(:,1:fn);
end




end
