%% resave data
clear all
%close hidden all
%addpath(genpath('D:\Box\Box\Grants\Cardiac DCE for 7 mins viability scan\Data\code\image-registration-master'));
%setup_image_registration_repository()
%% load Dicom
%subjectfolder='D:\Box\Box\Grants\Cardiac DCE for 7 mins viability scan\Data\code\RawData\Human\NCKU_May6\8 min T1 map\DCE\Slice3';
Rootfolder='C:\Users\yanghj\Box\Yang Lab\Projects\Cardiac Delayed Phase DCE\Data\Structured data\Animals';
Rootdir=dir(Rootfolder);
for r=[4]%[1 3 4 6 7 9 10];%1:(length(Rootdir)-2)
    filefolder=[Rootfolder,'\',Rootdir(r+2).name];
    filedir=dir(filefolder);

    for n=1:(length(filedir)-2)
        subjectfolder_in=[filefolder,'\',filedir(n+2).name];
        
        %DCE folder
        DCEfolder=[subjectfolder_in,'\DCE'];
        DCEdir=dir(fullfile(DCEfolder,'*.mat'))
        SeriesNumUpdated=0;
       % Run_Process_cardiacDCE_batch_func(subjectfolder_in)  
       for D=1:size(DCEdir)
           DCEMatpath=[DCEfolder,'\',DCEdir(D).name]
           load(DCEMatpath);
           SeriesNum=SeriesNumUpdated;
           Save_Files_DCE_batch;
           SeriesNumUpdated=SeriesNum;
       end
    end
end