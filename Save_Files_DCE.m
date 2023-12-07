
timepoints=uint16(ceil(time(2:end)/60))
imagedir=[subjectfolder,'\DCEnote\T',num2str(timepoints),'min\Images'];
dicominfo=DataPost(1).info;
%%
if SeriesNum<10
SeriesNum=savedicom(ECV,imagedir,dicominfo,SeriesNum,['ECV'])
SeriesNum=savedicom(dR1*1000,imagedir,dicominfo,SeriesNum,['dR1'])
end

%%
mkdir(imagedir);
SeriesNum=savedicom(F_cxm,imagedir,dicominfo,SeriesNum,['cmx_F_T',num2str(timepoints),'min'])
SeriesNum=savedicom(vp_cxm,imagedir,dicominfo,SeriesNum,['cxm_Vp_T',num2str(timepoints),'min'])
SeriesNum=savedicom(PS_cxm,imagedir,dicominfo,SeriesNum,['cxm_PS_T',num2str(timepoints),'min'])
SeriesNum=savedicom(ve_cxm,imagedir,dicominfo,SeriesNum,['cxm_Ve_T',num2str(timepoints),'min'])
SeriesNum=savedicom(Rsq_cxm,imagedir,dicominfo,SeriesNum,['cxm_Rsq_T',num2str(timepoints),'min'])
SeriesNum=savedicom(simLateR1cxm*1000,imagedir,dicominfo,SeriesNum,['SimR1cxm_permin_T',num2str(timepoints),'min'])%in sec
%savedicom(idximcxm,imagedir,dicominfo,'Kmeancluster_cxm')

SeriesNum=savedicom(tek_Ktrans,imagedir,dicominfo,SeriesNum,['etk_Ktrans_T',num2str(timepoints),'min'])
SeriesNum=savedicom(tek_Kep,imagedir,dicominfo,SeriesNum,['etk_Kep_T',num2str(timepoints),'min'])
SeriesNum=savedicom(tek_Vp,imagedir,dicominfo,SeriesNum,['etk_Vp_T',num2str(timepoints),'min'])
SeriesNum=savedicom(tek_Ve,imagedir,dicominfo,SeriesNum,['etk_Ve_T',num2str(timepoints),'min'])
SeriesNum=savedicom(tek_Rsq,imagedir,dicominfo,SeriesNum,['etk_Rsq_T',num2str(timepoints),'min'])
SeriesNum=savedicom(simLateR1etk*1000,imagedir,dicominfo,SeriesNum,['SimR1etk_permin_T',num2str(timepoints),'min'])%in sec
%savedicom(idximetk,imagedir,dicominfo,'Kmeancluster_etk')



mkdir([subjectfolder,'\DCE']);
save([subjectfolder,'\DCE\workspace T',num2str(timepoints),'min']) 
disp(['workspace is saved as ',subjectfolder,'\DCE\workspace'])

function [SeriesNum] = savedicom(var,imdir,dicominfo,SeriesNum,Name);

dicominfo.ProtocolName=Name;
dicominfo.SeriesDescription=Name;
dicominfo.SeriesNumber=dicominfo.SeriesNumber+1000+SeriesNum;
dicominfo.SeriesInstanceUID=[dicominfo.SeriesInstanceUID(1:end-1), num2str(SeriesNum+1)];



mat2dicom_DCE(var,[imdir,'\',Name,'\'],true,dicominfo);
disp([Name,' is saved as ',imdir,'\',Name])
SeriesNum=SeriesNum+1;
end
