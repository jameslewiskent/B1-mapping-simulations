function [settings] = load_sequence_settings(settings)
% Load settings for a particular sequence
if strcmpi(settings.Scheme,'SatTFL')
    [settings] = sattfl_settings(settings);
    
elseif strcmpi(settings.Scheme,'Sandwich')
    [settings] = sandwich_settings(settings);
    
elseif strcmpi(settings.Scheme,'SA2RAGE')
    [settings] = sa2rage_settings(settings);
    
elseif strcmpi(settings.Scheme,'AFI')
    [settings] = afi_settings(settings);
    
elseif strcmpi(settings.Scheme,'DREAM')
    [settings] = dream_settings(settings);
    
elseif strcmpi(settings.Scheme,'GRE')
    [settings] = gre_settings(settings);
else
    error('ABORTED: Scheme not recognised, please input either ''SatTFL'', ''Sandwich'', ''DREAM'', ''AFI'', ''SA2RAGE'' OR ''ALL''.')
end
end

