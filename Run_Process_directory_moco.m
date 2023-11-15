clear all
close hidden all

Rootfolder = '/Users/mona/Library/CloudStorage/Box-Box/Mona/patient_registration_LiTing_Mona';

subject = '075';
slice = 'CE8';
subjectfolder_in = fullfile(Rootfolder, subject);

PostconmocoFoldername = sprintf("%s/registration", subjectfolder_in);
% PostconmocoFoldername=[subjectfolder_in,'/PostconT1'];
if exist(PostconmocoFoldername,'dir')    
    Run_Process_cardiacDCE_unit_Moco_func(subjectfolder_in, 0, slice)  
end
