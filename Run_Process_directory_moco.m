clear all
close hidden all
addpath(genpath('.'))
Rootfolder = '/Users/mona/Library/CloudStorage/Dropbox/0.MAC-SYNC/0.PROJECT/Registration/Data/Liting_patient/T1MAPTEST';
Basefolder = '/Users/mona/Library/CloudStorage/Dropbox/0.MAC-SYNC/0.PROJECT/Registration/Data/Liting_patient/T1 maps for test ORIGINAL';
subject = '031';
slice = 'CE6_physical_model';
raw = 'MULTILEVEL';
subjectfolder_in = fullfile(Rootfolder, subject);
basefolder_in = fullfile(Basefolder, subject);
PostconmocoFoldername = fullfile(subjectfolder_in, 'REGISTRATION', raw);
% PostconmocoFoldername=[subjectfolder_in,'/PostconT1'];
tips = 'TV';
if exist(PostconmocoFoldername,'dir')    
    Run_Process_cardiacDCE_unit_Moco_func_TV(subjectfolder_in, basefolder_in, raw, slice, subject, tips);
    % Run_Process_cardiacDCE_unit_Moco_func(subjectfolder_in, basefolder_in, raw, slice, subject) 
end
