clear all
close hidden all

Rootfolder = 'C:\\Users\xinli\Box\Mona\patient_registration_LiTing_Mona';

subject = '017';
slice = 'CE5';
subjectfolder_in = fullfile(Rootfolder, subject);

PostconmocoFoldername = fullfile(subjectfolder_in, "registration");
% PostconmocoFoldername=[subjectfolder_in,'/PostconT1'];
if exist(PostconmocoFoldername,'dir')    
    Run_Process_cardiacDCE_unit_Moco_func(subjectfolder_in, 0, slice)  
end
