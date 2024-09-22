%
%function [FpFmZ,RR] = epg_rf(FpFmZ,alpha,phi)
%	Propagate EPG states through an RF rotation of 
%	alpha, with phase phi (both radians).
%	
%	INPUT:
%		FpFmZ = 3xN vector of F+, F- and Z states.
%		alpha = flip angle in radians.
%		phi = angle of rotation axis from Mx (radians).
%
%       OUTPUT:
%       FpFmZ = Updated FpFmZ state.
%		RR = RF rotation matrix (3x3).
%
%	SEE ALSO:
%		epg_grad, epg_grelax
%
%	B.Hargreaves.
%
function [FpFmZ,RR,Mag_Track] = epg_rf(FpFmZ,alpha,phi,RFTime,Mag_Track,settings)
% -- From Weigel at al, JMR 205(2010)276-285, Eq. 8.

if size(FpFmZ,2) > 10
[FpFmZ] = epg_trim(FpFmZ,settings.EPG_trim_threshold);
end

if (abs(alpha)>2*pi)
    warning('epg_rf:  Flip angle should be in radians!'); 
end

if (nargin < 3) 
    warning('Rotation axis not specified - assuming -My');
	phi=-pi/2; 
end

% -- Rotation matrix
RR = [(cos(alpha/2))^2 exp(2*1i*phi)*(sin(alpha/2))^2 -1i*exp(1i*phi)*sin(alpha);
      exp(-2*1i*phi)*(sin(alpha/2))^2 (cos(alpha/2))^2 1i*exp(-1i*phi)*sin(alpha);
      -1i/2*exp(-1i*phi)*sin(alpha) 1i/2*exp(1i*phi)*sin(alpha)      cos(alpha)];


FpFmZ = RR * FpFmZ;

% JK added for tracking magnetisation
if any(settings.Mag_Track_Flags == 1)
N_Samples = ceil(RFTime./settings.Mag_Track_dt); % In the case where Mag_Track_dt > RFTime, minimum 1 sample
mz = sum(epg_FZ2mz(FpFmZ),2);
Mag_Track(:,end + (1:N_Samples)) = [ones(1,N_Samples)*mz;Mag_Track(2,end) + cumsum(ones(1,N_Samples)*(RFTime/N_Samples),2)];
end

end


