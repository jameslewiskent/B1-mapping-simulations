function [already_ran,filename] = check_if_already_ran(settings)
% Check if simulations have already been run with the given parameters
filename = settings.filename;

% remove filename from settings as this depends on when it is ran as it contains date and time
try
    settings = rmfield(settings,'filename');
end

try
    % Also remove loop field name and values because this is irrelevant
    settings = rmfield(settings,'LoopFieldName');
    settings = rmfield(settings,'LoopValues');
    settings = rmfield(settings,'Additional_Loop_Counter');
end

try
    % Also remove loop field name and values because this is irrelevant
    settings = rmfield(settings,'LoopFieldName2');
    settings = rmfield(settings,'LoopValues2');
    settings = rmfield(settings,'Additional_Loop_Counter2');
end

try
    % Also ignore magtrack flags and temporal resolution of mag tracking as
    % it doesn't affect simulation results
    settings = rmfield(settings,'Mag_Track_Flags');
    settings = rmfield(settings,'Mag_Track_dt');
end

list = dir(settings.filepath);
list = list(~ismember({list.name},{'.','..'})); % Remove dirs
list = list(contains({list.name},'Results')); % Remove non-Results

already_ran = false;
for n = 1:length(list)
    %data = quickLoad(fullfile(settings.filepath,list(n).name),'settings');
    data = load(fullfile(settings.filepath,list(n).name),'settings');
    
    % remove filename from settings as this depends on when it is ran as it contains date and time
    try
        data.settings = rmfield(data.settings,'filename');
    end
    
    try
        data.settings = rmfield(data.settings,'LoopFieldName');
        data.settings = rmfield(data.settings,'LoopValues');
        data.settings = rmfield(data.settings,'Additional_Loop_Counter');
    end
    
    try
        data.settings = rmfield(data.settings,'LoopFieldName2');
        data.settings = rmfield(data.settings,'LoopValues2');
        data.settings = rmfield(data.settings,'Additional_Loop_Counter2');
    end
    
    try
        % Also ignore magtrack flags and temporal resolution of mag tracking as
        % it doesn't affect simulation results
        data.settings = rmfield(data.settings,'Mag_Track_Flags');
        data.settings = rmfield(data.settings,'Mag_Track_dt');
    end
    
    if isequaln(settings,data.settings) % Is equal (treat nans as equal)
        filename = list(n).name;
        already_ran = true;
        break
    else
        % Uncomment for debugging
%         %Find which field names don't match
%                         list_fieldnames = fieldnames(settings);
%                         for m = 1:length(list_fieldnames)
%                             if ~iscell(settings.(char(list_fieldnames(m))))
%                                 if ~all(isequaln(settings.(char(list_fieldnames(m))),data.settings.(char(list_fieldnames(m)))),'all')
%                                     disp(['Field name which does not match is: "', char(list_fieldnames(m)),'".']);
%                                     %disp(['Field name which does not match is: "', char(list_fieldnames(m)),'". Values are: ',num2str(settings.(char(list_fieldnames(m)))),' and ',num2str(data.settings.(char(list_fieldnames(m))))]);
%                                 end
%                             end
%                         end
%                         disp(' ')
                        
    end
    
    
end
end

