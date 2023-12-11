function [feedbak] =mat2dicom_DCE(im_Nonbin,save_fname,Dicominfo, info)%
%input matrix and directory name
%dicom info is true/flase to include dicom info or not
mkdir(save_fname);
counter = size(im_Nonbin,3);
orig = double(im_Nonbin);
 im_Nonbin=im_Nonbin/65535*1e2;%Decide to the percision point
for n = 1:counter,  
    if(n<10)
        niftiwrite(squeeze(orig(:,:,n)),strcat(save_fname,'00',num2str(n),'.nii'))
        dicomwrite(squeeze(im_Nonbin(:,:,n)),strcat(save_fname,'00',num2str(n),'.ima'),'CompressionMode','None','Quality',100);
    elseif (n < 99)
        niftiwrite(squeeze(orig(:,:,n)),strcat(save_fname,'0',num2str(n),'.nii'))
        dicomwrite(squeeze(im_Nonbin(:,:,n)),strcat(save_fname,'0',num2str(n),'.ima'),'CompressionMode','None','Quality',100);
    else
        niftiwrite(squeeze(orig(:,:,n)),strcat(save_fname,num2str(n),'.nii'))
        dicomwrite(squeeze(im_Nonbin(:,:,n)),strcat(save_fname,num2str(n),'.ima'),'CompressionMode','None','Quality',100);
    end
    if Dicominfo
       if(n<10)
            targetfile=strcat(save_fname,'00',num2str(n),'.ima');
       elseif (n < 99)
           targetfile=strcat(save_fname,'0',num2str(n),'.ima');
       else
           targetfile=strcat(save_fname,num2str(n),'.ima');
       end


%             
%                 I = dicomread(targetfile);
%                    dicomwrite(I,targetfile,info)
        
    end
end
feedbak=true;