function [Mz] = epg_FZ2mz(FpFmZ)

Ns = size(FpFmZ,2);	% Number of EPG states (N should be 2*Q-1 or more)
N = 2*Ns-1; 	% Minimum value for N.
		
% -- Use a matrix for FFT to support arbitrary N.
%  	This is because we are going from a discrete set of 
%	coefficients to a "continuous" distribution of Mx,My,Mz
x = (0:N-1)/N-0.5;	% Fraction of a cycle for spins locations.
ph = exp(1i*2*pi*x.'*(0:Ns-1));		% Phasors for Mz transform
FpFmZ(3,1)=FpFmZ(3,1)/2;		% Account for discretization
Mz = 2*real(ph*FpFmZ(3,:).');		% Transform to Mz

Mz = Mz.'/N;	% Form output vector.
end


