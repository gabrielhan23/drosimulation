% Load the dicom file for Quantitative Perfusion Study
% We sample different timepoints from different timesteps

% Mona, April 15, 2024

roi_path = fullfile(outputfolder, 'MASK', [label, '_ROI.mat']);
if exist(roi_path,'file')
    load(roi_path)
end


[V,spatial,dim] = dicomreadVolume(fullfile(subjectfolder, label));


files = natsortfiles(dir(fullfile(subjectfolder, label, 'MR*')));

assert(length(files) == size(V, 4), 'Loaded dicom volume and files mismatch')

% Load the pre contrast image, take the 4th frame as an example
idx = 4;
DataPre(1).img = V(:, :, 1, idx); 
DataPre(1).info = dicominfo(fullfile(files(idx).folder, files(idx).name));

% Load the extra as post contrast image, sample points
timeinstr = DataPre(1).info.AcquisitionTime;
timeino_pre(1)=str2num(timeinstr(1:2))*60*60+str2num(timeinstr(3:4))*60+str2num(timeinstr(5:end));%s

% Load the post contrast image, take the 3rd frame as an example
for i = 1:size(V, 4)
    DataPost_moco(i).img = dicomread(fullfile(files(i).folder, files(i).name)); 
    DataPost_moco(i).info = dicominfo(fullfile(files(i).folder, files(i).name));
    
    % Load the extra as post contrast image, sample points
    timeinstr = DataPost_moco(i).info.AcquisitionTime;
    timeino_moco(i)=str2num(timeinstr(1:2))*60*60+str2num(timeinstr(3:4))*60+str2num(timeinstr(5:end)); %s
end

%% Load the AIF LR

[V,spatial,dim] = dicomreadVolume(fullfile(subjectfolder, lr_label));


files = natsortfiles(dir(fullfile(subjectfolder, lr_label, 'MR*')));

assert(length(files) == size(V, 4), 'Loaded dicom volume and files mismatch')

idx = 4;
DataPre_aif(1).img = V(:, :, 1, idx); 
DataPre_aif(1).info = dicominfo(fullfile(files(idx).folder, files(idx).name));
% Load the post contrast image, take the 3rd frame as an example
for i = 1:size(V, 4)
    DataPost_aif_moco(i).img = dicomread(fullfile(files(i).folder, files(i).name)); 
    DataPost_aif_moco(i).info = dicominfo(fullfile(files(i).folder, files(i).name));
    
    % Load the extra as post contrast image, sample points
    timeinstr = DataPost_aif_moco(i).info.AcquisitionTime;
    timeino_aif_moco(i)=str2num(timeinstr(1:2))*60*60+str2num(timeinstr(3:4))*60+str2num(timeinstr(5:end)); %s
end

% aif and HR whole heart, the beginning 398 frames has same acquisition
% time

PreT1 = double(DataPre.img);
PreT1_aif = double(DataPre_aif.img);
[PostT1Reg, timein] = sortbyAcquisitionTime(DataPost_moco(1:length(timeino_aif_moco)), timeino_moco(1:length(timeino_aif_moco)));
[PostT1Reg_aif, timein_aif] = sortbyAcquisitionTime(DataPost_aif_moco, timeino_aif_moco);


function [Postimgs, timein] = sortbyAcquisitionTime(data, time)

for n = 1:length(data)
    Postimgs(:,:,n) = double(data(n).img);
end

% use the frame 4 time as 0
time = time-time(4);
% Sort time and image
[timesort,To] = sort(time);
Postimgs = Postimgs(:,:,To);
timein = timesort(4:end);
Postimgs = Postimgs(:, :, 4:end);
end