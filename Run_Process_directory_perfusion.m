clear all
close hidden all
addpath(genpath('.'))
Rootfolder = '/Users/mona/Library/CloudStorage/Box-Box/Quantitative Perfusion study/Human data';
subject = 'MR_2024-04-09_3541-7954-2013-0143_raw';
label = 'S37';
subjectfolder = fullfile(Rootfolder, subject);

output = '/Users/mona/Library/CloudStorage/Dropbox/0.MAC-SYNC/0.PROJECT/DCE/Data/QuantitativePerfusionStudy';
% Run_Process_cardiacDCE_unit_perfusion_func(subjectfolder, output, subject, label);
Run_Process_cardiacDCE_unit_perfusion_func;
    