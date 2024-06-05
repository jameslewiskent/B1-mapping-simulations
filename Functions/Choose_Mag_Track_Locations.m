function [settings] = Choose_Mag_Track_Locations(settings)
% Pick locations to track magnetisation

nPoints = input('How many locations would you like to track? ');
pts = readPoints(settings.Synthetic_T1s,nPoints);
settings.Mag_Track_SynInd = fliplr(round(pts)');

end

