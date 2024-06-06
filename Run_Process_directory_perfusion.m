clear all
close hidden all
addpath(genpath('.'))
Rootfolder = '/Users/mona/Library/CloudStorage/Dropbox/0.MAC-SYNC/0.PROJECT/DCE/Data/QuantitativePerfusionStudy/PatientData_052424';
% Rootfolder = 'C:\Users\xinli\Dropbox\0.MAC-SYNC\0.PROJECT\DCE\Data\QuantitativePerfusionStudy\PatientData_052424';
subject = 'MR_2024-05-24_0741-7954-8413-0106_raw';
label = 'S182';
lr_label = 'S132';
subjectfolder = fullfile(Rootfolder, subject);

resume = false;
% output = 'C:\Users\xinli\Dropbox\0.MAC-SYNC\0.PROJECT\DCE\Data\QuantitativePerfusionStudy_results';
output = '/Users/mona/Library/CloudStorage/Dropbox/0.MAC-SYNC/0.PROJECT/DCE/Data/QuantitativePerfusionStudy_results';
% Run_Process_cardiacDCE_unit_perfusion_func(subjectfolder, output, subject, label);
Run_Process_cardiacDCE_unit_perfusion_func;
    