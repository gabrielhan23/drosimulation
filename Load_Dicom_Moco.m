%% load Dicom
% Mona: change the load of data using index.

roi_path = fullfile(subjectfolder, 'DCE', label, 'ROI');
if exist(roi_path,'dir')
    load(fullfile(roi_path, 'ROI.mat'))
end


PreconFoldername = fullfile(subjectfolder, 'PreconT1');
%%
% Mona: change the postcon folder
PostconFoldername = fullfile(subjectfolder, 'PostconT1');

listPost = natsortfiles(dir(sprintf('%s/*T1*%s*', PostconFoldername, label)));
for n=1:length(listPost)
    fprintf(sprintf("Load the post contrast files: %s\n", fullfile(listPost(n).folder, listPost(n).name)))
    DataPost(n)=loaddicom(fullfile(listPost(n).folder, listPost(n).name));
    timeinstr{n}=DataPost(n).info.AcquisitionTime;
    timeino(n)=str2num(timeinstr{n}(1:2))*60*60+str2num(timeinstr{n}(3:4))*60+str2num(timeinstr{n}(5:end));%s
end

%%
listPre = natsortfiles(dir(sprintf('%s/*T1*%s*', PreconFoldername, label)));
for n=1:length(listPre)
    fprintf(sprintf("Load the pre contrast files: %s\n", fullfile(listPre(n).folder, listPre(n).name)))
    DataPre(n)=loaddicom(fullfile(listPre(n).folder, listPre(n).name));
end
%%
% subjectfolder_in/registration/CE4/single_exp/stage2/T1_2MINPOSTCE4_MAP_T1_Mona
PostconmocoFoldername = fullfile(subjectfolder_in, 'registration', label, 'single_exp/stage2');
if exist(PostconmocoFoldername,'dir')
  
    listPost = natsortfiles(dir(sprintf('%s/*T1*', PostconmocoFoldername)));
    for n=1:length(listPost)
        fprintf(sprintf("Load the post contrast motion corrected files: %s\n", fullfile(listPost(n).folder, listPost(n).name)))
        DataPost_moco(n) = loaddicom(fullfile(listPost(n).folder, listPost(n).name));
        timeinstr_moco{n}=DataPost_moco(n).info.AcquisitionTime;
        timeino_moco(n)=str2num(timeinstr_moco{n}(1:2))*60*60+...
            str2num(timeinstr_moco{n}(3:4))*60+str2num(timeinstr_moco{n}(5:end));%s
    end
end
%%
function data = loaddicom(path)
% load the dicom files
list=dir(sprintf('%s/*.dcm', path));

% Mona: change due to exist of .DS_Store
data.img=dicomread(fullfile(list(end).folder, list(end).name));
data.info=dicominfo(fullfile(list(end).folder, list(end).name));

end