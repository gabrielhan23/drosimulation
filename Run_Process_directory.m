

clear all
close hidden all
%addpath(genpath('D:\Box\Box\Grants\Cardiac DCE for 7 mins viability scan\Data\code\image-registration-master'));
%setup_image_registration_repository()
%% load Dicom
%subjectfolder='D:\Box\Box\Grants\Cardiac DCE for 7 mins viability scan\Data\code\RawData\Human\NCKU_May6\8 min T1 map\DCE\Slice3';
Rootfolder='C:\Users\yanghj\Box\Yang Lab\Projects\Cardiac Delayed Phase DCE\Data\Structured data\Animals';
%Rootfolder='D:\Box\Box\Yang Lab\Projects\Cardiac Delayed Phase DCE\Data\Structured data\Animals';
Rootdir=dir(Rootfolder);
%%
for r=7:(length(Rootdir)-2)
    filefolder=[Rootfolder,'\',Rootdir(r+2).name];
    filedir=dir(filefolder);

    for n=1:(length(filedir)-2)
        subjectfolder_in=[filefolder,'\',filedir(n+2).name]
        %Run_Process_cardiacDCE_batch_func(subjectfolder_in)  
        Run_Process_cardiacDCE_unit_func(subjectfolder_in)  
    end
end

%next[11]