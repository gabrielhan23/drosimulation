
timepoints=uint16(ceil(time(2:end)/60))
imagedir=[subjectfolder,'\DCE\T',num2str(timepoints),'min\Images'];
dicominfo=DataPost(1){ns}.info{ns};

mkdir(imagedir);
savedicom(F,imagedir,dicominfo,'cmx_F')
savedicom(vp,imagedir,dicominfo,'cxm_Vp')
savedicom(PS,imagedir,dicominfo,'cxm_PS')
savedicom(ve,imagedir,dicominfo,'cxm_Ve')
savedicom(fitresultsDCEcxm(:,:,5),imagedir,dicominfo,'cxm_Rsq')
savedicom(simLateR1cxm*1000,imagedir,dicominfo,'SimR1cxm_permin')%in sec
%savedicom(idximcxm,imagedir,dicominfo,'Kmeancluster_cxm')

savedicom(tek_Ktrans,imagedir,dicominfo,'etk_Ktrans')
savedicom(tek_Kep,imagedir,dicominfo,'etk_Kep')
savedicom(tek_Ve,imagedir,dicominfo,'etk_Ve')
savedicom(tek_Rsq,imagedir,dicominfo,'etk_Rsq')
savedicom(simLateR1etk*1000,imagedir,dicominfo,'SimR1etk_permin')%in sec
%savedicom(idximetk,imagedir,dicominfo,'Kmeancluster_etk')

savedicom(ECV,imagedir,dicominfo,'ECV')

mkdir([subjectfolder,'\DCE']);
save([subjectfolder,'\DCE\T',num2str(timepoints),'min\workspace']) 
disp(['workspace is saved as ',subjectfolder,'\DCE\T',num2str(timepoints),'min\workspace'])

function [] = savedicom(var,imdir,dicominfo,Name);

dicominfo.ProtocolName=Name;
mat2dicom_DCE(var,[imdir,'\',Name,'\'],true,dicominfo);
disp([Name,' is saved as ',imdir,'\',Name])
end
