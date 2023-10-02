function [settings] = load_sequence_settings(settings)
% Load settings for a particular sequence
if strcmp(settings.Scheme,'SatTFL')
    [settings] = sattfl_settings(settings);
    
elseif strcmp(settings.Scheme,'Sandwich')
    [settings] = sandwich_settings(settings);
    
elseif strcmp(settings.Scheme,'SA2RAGE')
    [settings] = sa2rage_settings(settings);
    
elseif strcmp(settings.Scheme,'AFI')
    [settings] = afi_settings(settings);
    
elseif strcmp(settings.Scheme,'DREAM')
    [settings] = dream_settings(settings);
else
    error('ABORTED: Scheme not recognised, please input either ''SatTFL'', ''Sandwich'', ''DREAM'', ''AFI'', ''SA2RAGE'' OR ''ALL''.')
end
end

