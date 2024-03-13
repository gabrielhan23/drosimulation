clear all
close hidden all
addpath(genpath('.'))
Rootfolder = 'C:\Users\xinli\Dropbox\0.MAC-SYNC\0.PROJECT\Registration\Data\Liting_patient\T1MAPTEST';
Basefolder = 'C:\Users\xinli\Dropbox\0.MAC-SYNC\0.PROJECT\Registration\Data\Liting_patient\T1 maps for test ORIGINAL';
subject = '031';
slice = 'CE6';
raw = 'MULTILEVEL';
subjectfolder_in = fullfile(Rootfolder, subject);
basefolder_in = fullfile(Basefolder, subject);
PostconmocoFoldername = fullfile(subjectfolder_in, 'REGISTRATION', raw);
% PostconmocoFoldername=[subjectfolder_in,'/PostconT1'];
if exist(PostconmocoFoldername,'dir')    
    Run_Process_cardiacDCE_unit_Moco_func(subjectfolder_in, basefolder_in, raw, slice, subject)  
end
