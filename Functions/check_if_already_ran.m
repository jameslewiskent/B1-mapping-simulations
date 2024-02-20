function [already_ran,savefilename,filename] = check_if_already_ran(settings)
% Check if simulations have already been run with the given parameters
savefilename = settings.savefilename;

% remove filename from settings as this depends on when it is ran as it contains date and time
try
settings = rmfield(settings,'savefilename'); settings = rmfield(settings,'filename');
end

list = dir(settings.filepath);
list = list(~ismember({list.name},{'.','..'})); % Remove dirs
list = list(contains({list.name},'Results')); % Remove non-Results

already_ran = false;
for n = 1:length(list)
    data = load(fullfile(settings.filepath,list(n).name),'settings');
    
    % remove filename from settings as this depends on when it is ran as it contains date and time
    try
    data.settings = rmfield(data.settings,'savefilename'); 
    filename = data.settings.filename;
    data.settings = rmfield(data.settings,'filename');
    end
    
    if isequaln(settings,data.settings) % Is equal (treat nans as equal)
        savefilename = list(n).name;
        already_ran = true;
        break
    end
end
end

