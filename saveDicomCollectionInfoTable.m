% Get the Dicom Collection Information

DICOM_FOLDER = 'C:\Users\xinli\Dropbox\0.MAC-SYNC\0.PROJECT\DCE\Data\QuantitativePerfusionStudy\PatientData_052424\MR_2024-05-24_0741-7954-8413-0106_raw';
details = dicomCollection(DICOM_FOLDER);

% save the details in xlsx sheet
writetable(details,fullfile(DICOM_FOLDER, 'summary_table.csv'), 'Delimiter',',','QuoteStrings','all', 'WriteRowNames', true);

%% Load the dicom volume
[V,spatial,dim] = dicomreadVolume(fullfile(DICOM_FOLDER, 'S182'));


files = natsortfiles(dir(fullfile(DICOM_FOLDER, 'S182/MR*')));

%% Sample the images

% Load the pre contrast image, take the 3rd frame as an example
idx = 3;
DataPre(1).data = V(:, :, 1, idx); 
DataPre(1).info = dicominfo(fullfile(files(idx).folder, files(idx).name));

% Load the extra as post contrast image, sample points
timeinstr = DataPre(1).info.AcquisitionTime;
timeino_pre(1)=str2num(timeinstr(1:2))*60*60+str2num(timeinstr(3:4))*60+str2num(timeinstr(5:end));%s

% Load the post contrast image, take the 3rd frame as an example
for i = 1:size(V, 4)
    DataPost_moco(i).data = V(:, :, 1, i); 
    DataPost_moco(i).info = dicominfo(fullfile(files(i).folder, files(i).name));
    
    % Load the extra as post contrast image, sample points
    timeinstr = DataPre(i).info.AcquisitionTime;
    timeino_moco(i)=str2num(timeinstr(1:2))*60*60+str2num(timeinstr(3:4))*60+str2num(timeinstr(5:end)); %s
end