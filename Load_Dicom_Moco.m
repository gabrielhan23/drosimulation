%% load Dicom
% Mona: change the load of data using index.

% add the aciquisition time
load(fullfile(subjectfolder, 'acquisitionTime.mat'))

roi_path = fullfile(subjectfolder, 'MASK', [label(1:3), '_ROI.mat']);
if exist(roi_path,'file')
    load(roi_path)
end


PreconFoldername = fullfile(subjectfolder, 'REGISTRATION', raw, label, 'PRECONT1');

% Mona: change the postcon folder
PostconFoldername = fullfile(subjectfolder, 'REGISTRATION', raw, label, 'POSTCONT1');

listPost = natsortfiles(dir(fullfile(PostconFoldername, '*T1MAP_MONA')));
for n=1:length(listPost)
    disp(sprintf("Load the post contrast files: %s", fullfile(listPost(n).folder, listPost(n).name)))

    tmp = regexp(listPost(n).name, '^(.*_MAP)', 'match');
    name = tmp{1};
    
    DataPost_moco(n) = loaddicom(fullfile(listPost(n).folder, listPost(n).name), basefolder, extractAfter(name, [subject, '_']));
    timeino_moco(n) = timeino_map(name);
    % 
    % % name = listPost(n).name;
    % % info_folder = fullfile(subjectfolder, extractAfter(extractBefore(name, '_T1MAP'), [subject, '_']));
    % % orig_dicoms = dir(info_folder);
    % % info = dicominfo(fullfile(orig_dicoms(end).folder, orig_dicoms(end).name));
    % DataPost(n) = loaddicom(fullfile(listPost(n).folder, listPost(n).name));
    % timeinstr{n}=DataPost(n).info.AcquisitionTime;
    % timeino(n)=str2num(timeinstr{n}(1:2))*60*60+str2num(timeinstr{n}(3:4))*60+str2num(timeinstr{n}(5:end));%s
end


listPre = natsortfiles(dir(fullfile(PreconFoldername, '*T1MAP_MONA')));
for n=1:length(listPre)
    disp(sprintf("Load the post contrast files: %s", fullfile(listPre(n).folder, listPre(n).name)))

    tmp = regexp(listPre(n).name, '^(.*_MAP)', 'match');
    name = tmp{1};
    
    DataPre(n)=loaddicom(fullfile(listPre(n).folder, listPre(n).name), basefolder, extractAfter(name, [subject, '_']));
    timeino_pre(n) = timeino_map(name);
end

%%
% subjectfolder_in/registration/CE4/single_exp/stage2/T1_2MINPOSTCE4_MAP_T1_Mona
% if raw == 1
%     PostconmocoFoldername = fullfile(subjectfolder_in, 'registration', label, 'single_exp', 'stage1');
% else
%     PostconmocoFoldername = fullfile(subjectfolder_in, 'registration', label, 'single_exp', 'stage2');
% end
% PostconmocoFoldername = fullfile()
% if exist(PostconmocoFoldername,'dir')
% 
%     listPost = natsortfiles(dir(fullfile(PostconmocoFoldername, '*POST*')));
%     for n=1:length(listPost)
%         disp(sprintf("Load the post contrast motion corrected files: %s", fullfile(listPost(n).folder, listPost(n).name)))
%         DataPost_moco(n) = loaddicom(fullfile(listPost(n).folder, listPost(n).name));
%         timeinstr_moco{n}=DataPost_moco(n).info.AcquisitionTime;
%         timeino_moco(n)=str2num(timeinstr_moco{n}(1:2))*60*60+...
%             str2num(timeinstr_moco{n}(3:4))*60+str2num(timeinstr_moco{n}(5:end));%s
%     end
%     timeino_moco = timeino;
% end

%%
function data = loaddicom(path, basefolder, name)
% load the dicom files
list=dir(fullfile(path, '*'));

% Mona: change due to exist of .DS_Store
data.img=niftiread(fullfile(list(end).folder, list(end).name));
% data.img = dicomread(fullfile(list(end).folder, list(end).name));
list=dir(fullfile(basefolder, name, '*'));
if ~isempty(list)
    data.info=dicominfo(fullfile(list(end).folder, list(end).name));
else
    data.info = [];
end
end