function settings = Create_Filename(settings)
% Generate filename

settings.filename = ['Results_',char(datetime('now','TimeZone','local','Format','d-MMM-y-HH-mm'))];

n = 1;
alf = 'a':'z';
while isfile(fullfile(settings.filepath,settings.filename))
    % File already exists, potentially ran in same minute, give suffix
    settings.filename = ['Results_',char(datetime('now','TimeZone','local','Format','d-MMM-y-HH-mm')),alf(n)];
    n = n + 1;
end

end

