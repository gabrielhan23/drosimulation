%% load Dicom
%subjectfolder='D:\Box\Box\Grants\Cardiac DCE for 7 mins viability scan\Data\code\RawData\Lisbon_D7';
if exist([subjectfolder,'\DCE\workspace.mat'])
    load([subjectfolder,'\DCE\workspace.mat'])
else
% subjectfolder='D:\Box\Box\Grants\Cardiac DCE for 7 mins viability scan\Data\code\Data\SofiaD7';
PreconFoldername=[subjectfolder,'\PreconT1'];
PostconFoldername=[subjectfolder,'\PostconT1'];

listPost=dir(PostconFoldername);
for n=1:length(listPost)-2
    DataPost(n)=loaddicom([PostconFoldername,filesep,listPost(n+2).name]);
    timeinos{n}=DataPost(n).timeino;
end
listPre=dir(PreconFoldername);
for n=1:length(listPre)-2
    DataPre(n)=loaddicom([PreconFoldername,filesep,listPre(n+2).name]);
end
end

%%
function data = loaddicom(path)
list=dir(path);
%if ~(strcmp('IMA',list(3).name(end-2:end)))

%  fprintf('no Dicom in the directory')  
%else
for c=1:length(list)-2
data.img{c}=dicomread([path,filesep,list(c+2).name]);
data.info{c}=dicominfo([path,filesep,list(c+2).name]);
data.timeinstr{c}=data.info{c}.AcquisitionTime;
data.timeino(c)=str2num(data.timeinstr{c}(1:2))*60*60+str2num(data.timeinstr{c}(3:4))*60+str2num(data.timeinstr{c}(5:end));%s

end

%end
end