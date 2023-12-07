%% load Dicom
%subjectfolder='D:\Box\Box\Grants\Cardiac DCE for 7 mins viability scan\Data\code\RawData\Lisbon_D7';
DCEfolder=[subjectfolder,'\DCE'];
matfilename=dir([subjectfolder,'\DCE\*.mat']);
%if exist([subjectfolder,'\DCE\workspace.mat'])
if ~isempty(matfilename)
MannualMaskexist=0;
load([subjectfolder,'\DCE\',matfilename(1).name],'Blood', 'HeartMask', 'LVendo', 'Myomask', 'MIMask', 'MVOMask', 'RemoteMask')
if exist([subjectfolder,'\DCE\ROI_MannualMaskbig.mat'], 'file') 
load([subjectfolder,'\DCE\ROI_MannualMaskbig.mat']);
MannualMaskexist=1;
end
rmdir([subjectfolder,'\DCE\'],'s');
pause(1)
mkdir([subjectfolder,'\DCE\'],'s');

if MannualMaskexist
    save([DCEfolder,'\ROI_MannualMaskbig.mat'],'MannualMask');
end
    save([DCEfolder,'\ROI.mat'],'Blood', 'HeartMask', 'LVendo', 'Myomask', 'MIMask', 'MVOMask', 'RemoteMask');
end
% subjectfolder='D:\Box\Box\Grants\Cardiac DCE for 7 mins viability scan\Data\code\Data\SofiaD7';
PreconFoldername=[subjectfolder,'\PreconT1'];
PostconFoldername=[subjectfolder,'\PostconT1'];

listPost=dir(PostconFoldername);
for n=1:length(listPost)-2
    [PostconFoldername,filesep,listPost(n+2).name]
    DataPost(n)=loaddicom([PostconFoldername,filesep,listPost(n+2).name]);
    timeinstr{n}=DataPost(n).info.AcquisitionTime;
    timeino(n)=str2num(timeinstr{n}(1:2))*60*60+str2num(timeinstr{n}(3:4))*60+str2num(timeinstr{n}(5:end));%s
end
listPre=dir(PreconFoldername);
for n=1:length(listPre)-2
    DataPre(n)=loaddicom([PreconFoldername,filesep,listPre(n+2).name]);
end
%% load moco DICOM

PostconmocoFoldername=[subjectfolder,'\registration\stage2'];
if exist(PostconmocoFoldername,'dir')
listPost=dir(PostconmocoFoldername)
for n=1:length(listPost)-2
    [PostconmocoFoldername,filesep,listPost(n+2).name]
      DataPost_moco(n)=loaddicom([PostconmocoFoldername,filesep,listPost(n+2).name])
    timeinstr_moco{n}=DataPost_moco(n).info.AcquisitionTime;
    timeino_moco(n)=str2num(timeinstr_moco{n}(1:2))*60*60+...
        str2num(timeinstr_moco{n}(3:4))*60+str2num(timeinstr_moco{n}(5:end));%s
end
end

%%
function data = loaddicom(path)
list=dir(path);
%if ~(strcmp('IMA',list(3).name(end-2:end)))

%  fprintf('no Dicom in the directory')  
%else
data.img=dicomread([path,filesep,list(3).name]);
data.info=dicominfo([path,filesep,list(3).name]);

%end
end