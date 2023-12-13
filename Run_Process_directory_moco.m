clear all
close hidden all
addpath(genpath('.'))
Rootfolder = '/Users/mona/Library/CloudStorage/GoogleDrive-xinqili16@g.ucla.edu/My Drive/Data/Registration/patient_registration_LiTing_Mona';

subject = '075';
slice = 'CE4';
raw = 0;
subjectfolder_in = fullfile(Rootfolder, subject);

PostconmocoFoldername = fullfile(subjectfolder_in, "registration");
% PostconmocoFoldername=[subjectfolder_in,'/PostconT1'];
if exist(PostconmocoFoldername,'dir')    
    Run_Process_cardiacDCE_unit_Moco_func(subjectfolder_in, raw, slice)  
end
