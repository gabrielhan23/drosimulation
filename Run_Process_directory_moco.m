clear all
close hidden all
addpath(genpath('.'))
Rootfolder = 'C:\Users\xinli\Dropbox\Projects\rPCA-GroupRegNet\test_data\patient_registration_LiTing_Mona';

subject = '031';
slice = 'CE6';
raw = 1;
subjectfolder_in = fullfile(Rootfolder, subject);

PostconmocoFoldername = fullfile(subjectfolder_in, "registration");
% PostconmocoFoldername=[subjectfolder_in,'/PostconT1'];
if exist(PostconmocoFoldername,'dir')    
    Run_Process_cardiacDCE_unit_Moco_func(subjectfolder_in, raw, slice, subject)  
end
