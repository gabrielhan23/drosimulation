%% load Dicom
% Mona: change the load of data using index.

% add the aciquisition time
load(fullfile(subjectfolder, 'acquisitionTime.mat'))

roi_path = fullfile(subjectfolder, 'DCE', label, 'ROI');
if exist(roi_path,'dir')
    load(fullfile(roi_path, 'ROI.mat'))
end


PreconFoldername = fullfile(subjectfolder, 'PreconT1');

% Mona: change the postcon folder
PostconFoldername = fullfile(subjectfolder, 'PostconT1');

listPost = natsortfiles(dir(fullfile(PostconFoldername, sprintf('*%s*', label))));
for n=1:length(listPost)
    disp(sprintf("Load the post contrast files: %s", fullfile(listPost(n).folder, listPost(n).name)))
    DataPost(n)=loaddicom(fullfile(listPost(n).folder, listPost(n).name));

    
    name = listPost(n).name;
    info_folder = fullfile(subjectfolder, extractAfter(extractBefore(name, '_T1MAP'), [subject, '_']));
    orig_dicoms = dir(info_folder);
    info = dicominfo(fullfile(orig_dicoms(end).folder, orig_dicoms(end).name));
    DataPost(n).info = info;
%     timeinstr{n}=DataPost(n).info.AcquisitionTime;
%     timeino(n)=str2num(timeinstr{n}(1:2))*60*60+str2num(timeinstr{n}(3:4))*60+str2num(timeinstr{n}(5:end));%s
    timeino(n) = timeino_map(listPost(n).name);
end


%%
% subjectfolder_in/registration/CE4/single_exp/stage2/T1_2MINPOSTCE4_MAP_T1_Mona
if raw == 1
    PostconmocoFoldername = fullfile(subjectfolder_in, 'registration', label, 'single_exp', 'stage1');
else
    PostconmocoFoldername = fullfile(subjectfolder_in, 'registration', label, 'single_exp', 'stage2');
end
if exist(PostconmocoFoldername,'dir')
  
    listPost = natsortfiles(dir(fullfile(PostconmocoFoldername, '*POST*')));
    for n=1:length(listPost)
        disp(sprintf("Load the post contrast motion corrected files: %s", fullfile(listPost(n).folder, listPost(n).name)))
        DataPost_moco(n) = loaddicom(fullfile(listPost(n).folder, listPost(n).name));
%         timeinstr_moco{n}=DataPost_moco(n).info.AcquisitionTime;
%         timeino_moco(n)=str2num(timeinstr_moco{n}(1:2))*60*60+...
%             str2num(timeinstr_moco{n}(3:4))*60+str2num(timeinstr_moco{n}(5:end));%s
    end
    timeino_moco = timeino;
end

listPre = natsortfiles(dir(fullfile(PostconmocoFoldername, '*PRE*')));
for n=1:length(listPre)
    disp(sprintf("Load the pre contrast files: %s", fullfile(listPre(n).folder, listPre(n).name)))
    DataPre(n)=loaddicom(fullfile(listPre(n).folder, listPre(n).name));
end

%%
function data = loaddicom(path)
% load the dicom files
list=dir(fullfile(path, '*T1*.nii'));

% Mona: change due to exist of .DS_Store
data.img=niftiread(fullfile(list(end).folder, list(end).name));
data.img = data.img';
list=dir(fullfile(path, '*T1*.dcm'));
if ~isempty(list)
    data.info=dicominfo(fullfile(list(end).folder, list(end).name));
else
    data.info = [];
end
end