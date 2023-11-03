function [already_ran,filename] = check_if_already_ran(settings)
% Check if simulations have already been run with the given parameters

list = dir(settings.filepath);
list = list(~ismember({list.name},{'.','..'}));
list = list(~ismember({list.name},{settings.lookup_filename})); % Remove any lookup tables from the list

already_ran = false;
filename = 'NA';
for n = 1:length(list)
    data = load(fullfile(settings.filepath,list(n).name),'settings');
    previous_settings = data.settings;

    if isequal(settings,previous_settings)
        filename = list(n).name;
        already_ran = true;
        break
    end
end
end

