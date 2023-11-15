
timepoints=uint16(ceil(time(2:end)/60))
imagedir=[subjectfolder filesep folder, '_MOCO' filesep 'T', num2str(timepoints),'min' filesep 'Images'];


dicominfo=DataPost(1).info;

mkdir(imagedir);

%%
if SeriesNum<10
SeriesNum=savedicom(ECV,imagedir,dicominfo,SeriesNum,['ECV'])
temp=dR1(:,:,1);
SeriesNum=savedicom(dR1/median(temp(Blood))*100,imagedir,dicominfo,SeriesNum,['dR1'])
end
%%

%%
% scale Variable for Dicom range


vp_cxm_w=DerivedDCEMaps.vp_cxm;
F_cxm_w=DerivedDCEMaps.F_cxm*100;
PS_cxm_w=DerivedDCEMaps.PS_cxm*100;
ve_cxm_w=DerivedDCEMaps.ve_cxm;
simLateR1cxm_w=DerivedDCEMaps.simLateR1cxm*10;




vp_cxmD_w=DerivedDCEMaps.vp_cxmD;
PS_cxmD_w=DerivedDCEMaps.PS_cxmD*100;
ve_cxmD_w=DerivedDCEMaps.ve_cxmD;
simLateR1cxmD_w=DerivedDCEMaps.simLateR1cxmD*10;

tek_Ktrans_w=DerivedDCEMaps.tek_Ktrans*100;
tek_Kep_w=DerivedDCEMaps.tek_Kep*100;
tek_Vp_w=DerivedDCEMaps.tek_Vp;
tek_Ve_w=DerivedDCEMaps.tek_Ve;
simLateR1etk_w=DerivedDCEMaps.simLateR1etk*10;

%
vp_cxmmix_w=DerivedDCEMaps.vp_cxmmix;
F_cxmmix_w=DerivedDCEMaps.F_cxmmix*100;
PS_cxmmix_w=DerivedDCEMaps.PS_cxmmix*100;
ve_cxmmix_w=DerivedDCEMaps.ve_cxmmix;
simLateR1cxmmix_w=DerivedDCEMaps.simLateR1cxmmix*10;


%}
%%
SeriesNum=savedicom(F_cxm_w,imagedir,dicominfo,SeriesNum,['cmx_F_T',num2str(timepoints),'min'])
SeriesNum=savedicom(vp_cxm_w,imagedir,dicominfo,SeriesNum,['cxm_Vp_T',num2str(timepoints),'min'])
SeriesNum=savedicom(PS_cxm_w,imagedir,dicominfo,SeriesNum,['cxm_PS_T',num2str(timepoints),'min'])
SeriesNum=savedicom(ve_cxm_w,imagedir,dicominfo,SeriesNum,['cxm_Ve_T',num2str(timepoints),'min'])
SeriesNum=savedicom(Rsq_cxm,imagedir,dicominfo,SeriesNum,['cxm_Rsq_T',num2str(timepoints),'min'])
SeriesNum=savedicom(simLateR1cxm_w,imagedir,dicominfo,SeriesNum,['SimR1cxm_T',num2str(timepoints),'min'])%in sec
%savedicom(idximcxm,imagedir,dicominfo,'Kmeancluster_cxm')

SeriesNum=savedicom(vp_cxmD_w,imagedir,dicominfo,SeriesNum,['cxmD_Vp_T',num2str(timepoints),'min'])
SeriesNum=savedicom(PS_cxmD_w,imagedir,dicominfo,SeriesNum,['cxmD_PS_T',num2str(timepoints),'min'])
SeriesNum=savedicom(ve_cxmD_w,imagedir,dicominfo,SeriesNum,['cxmD_Ve_T',num2str(timepoints),'min'])
SeriesNum=savedicom(Rsq_cxmD,imagedir,dicominfo,SeriesNum,['cxmD_Rsq_T',num2str(timepoints),'min'])
SeriesNum=savedicom(simLateR1cxmD_w,imagedir,dicominfo,SeriesNum,['SimR1cxmD_T',num2str(timepoints),'min'])%in sec


SeriesNum=savedicom(tek_Ktrans_w,imagedir,dicominfo,SeriesNum,['etk_Ktrans_T',num2str(timepoints),'min'])
SeriesNum=savedicom(tek_Kep_w,imagedir,dicominfo,SeriesNum,['etk_Kep_T',num2str(timepoints),'min'])
SeriesNum=savedicom(tek_Vp_w,imagedir,dicominfo,SeriesNum,['etk_Vp_T',num2str(timepoints),'min'])
SeriesNum=savedicom(tek_Ve_w,imagedir,dicominfo,SeriesNum,['etk_Ve_T',num2str(timepoints),'min'])
SeriesNum=savedicom(tek_Rsq,imagedir,dicominfo,SeriesNum,['etk_Rsq_T',num2str(timepoints),'min'])
SeriesNum=savedicom(simLateR1etk_w,imagedir,dicominfo,SeriesNum,['SimR1etk_T',num2str(timepoints),'min'])%in sec
%savedicom(idximetk,imagedir,dicominfo,'Kmeancluster_etk')

%

SeriesNum=savedicom(F_cxmmix_w,imagedir,dicominfo,SeriesNum,['cmxmix_F_T',num2str(timepoints),'min'])
SeriesNum=savedicom(vp_cxmmix_w,imagedir,dicominfo,SeriesNum,['cxmmix_Vp_T',num2str(timepoints),'min'])
SeriesNum=savedicom(PS_cxmmix_w,imagedir,dicominfo,SeriesNum,['cxmmix_PS_T',num2str(timepoints),'min'])
SeriesNum=savedicom(ve_cxmmix_w,imagedir,dicominfo,SeriesNum,['cxmmix_Ve_T',num2str(timepoints),'min'])
SeriesNum=savedicom(Rsq_cxmmix,imagedir,dicominfo,SeriesNum,['cxmmix_Rsq_T',num2str(timepoints),'min'])
SeriesNum=savedicom(simLateR1cxmmix_w,imagedir,dicominfo,SeriesNum,['SimR1cxmmix_T',num2str(timepoints),'min'])%in sec


%% Modified for masked cxm
simLateR1cxm_masked_w=DerivedDCEMaps.simLateR1cxm_Masked*20;
SeriesNum=savedicom(simLateR1cxm_masked_w,imagedir,dicominfo,SeriesNum,['SimR1cxm_Masked_T',num2str(timepoints),'min'])%in sec
BloodMask_cxm=DerivedDCEMaps.BloodMask_cxm;
SeriesNum=savedicom(BloodMask_cxm,imagedir,dicominfo,SeriesNum,['BloodMask_cxm'])


%%

save(fullfile(subjectfolder, 'DCE', label, ['workspace T',num2str(timepoints),'min']))
disp(['workspace is saved as ', fullfile(subjectfolder, 'DCE', label, ['workspace T',num2str(timepoints),'min'])])



function [SeriesNum] = savedicom(var,imdir,dicominfo,SeriesNum,Name);
dicominfo=setdicominfo(var,dicominfo,SeriesNum,Name);
switch inputname(1)
     case {'simLateR1cxm_w','simLateR1etk_w','simLateR1cxm_masked_w', 'simLateR1cxmmix_w'}
       DicomTime=dicominfo.AcquisitionTime;
       for t=1:size(var,3)
        tempT=DicomTime; 
        tempmin=str2num(tempT(3:4))+t;
        temphr=str2num(tempT(1:2));

            temphr=temphr+floor(tempmin/60);
            tempmin=mod(tempmin,60);
            
        tempT(1:2)=num2str(temphr);
        tempT(3:4)=num2str(tempmin);

        Dicomtimearray{t}=tempT;
       end
       Dicomtimearray
    mat2dicom_DCE_sim(var,[imdir filesep Name filesep],true,dicominfo,Dicomtimearray);
    disp([Name,' is saved as ',imdir filesep Name])
    SeriesNum=SeriesNum+1; 
    otherwise
    mat2dicom_DCE(var,[imdir filesep Name filesep],true,dicominfo);
    disp([Name,' is saved as ',imdir filesep,Name])
    SeriesNum=SeriesNum+1;
end


end

function [dicominfo] = setdicominfo(var,dicominfo,SeriesNum,Name);

switch inputname(1)
    case 'F_cxm_w'
        dicominfo.SeriesDescription='Flow(L/100s/ml)';
    case 'vp_cxm_w'
        dicominfo.SeriesDescription='Blood volume percentage(1/10000)';
    case 'PS_cxm_w'
        dicominfo.SeriesDescription='CA exchange rate(1/100s)';
    case 've_cxm_w'
        dicominfo.SeriesDescription='EES volume percentage(1/10000)';
    case 'tek_Kep_w'
        dicominfo.SeriesDescription='Washout rate(1/100s)';
    case  'tek_Ktrans_w'
        dicominfo.SeriesDescription='Wash in rate(1/100s)';
    case 'tek_Ve_w'
        dicominfo.SeriesDescription='Washout ratio(Ktrans/Kep)(1/10000)';
    case 'tek_Vp_w'
        dicominfo.SeriesDescription='Blood volume percentage(1/10000)';
    case {'simLateR1etk_w', 'simLateR1cxm_w','simLateR1cxm_masked_w', 'simLateR1cxmmix_w'}
        dicominfo.SeriesDescription='sim per minute(mmol.s/l*10)';
        
       
    otherwise 
        dicominfo.SeriesDescription=Name;
end
    
    
dicominfo.ProtocolName=Name;

dicominfo.SeriesNumber=dicominfo.SeriesNumber+1000+SeriesNum;
dicominfo.SeriesInstanceUID=[dicominfo.SeriesInstanceUID(1:end-1), num2str(SeriesNum+1)];

end