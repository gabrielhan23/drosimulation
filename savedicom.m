function [SeriesNum] = savedicom(var,imdir,dicominfo,SeriesNum,Name)

dicominfo.ProtocolName=Name;
dicominfo.SeriesDescription=Name;
dicominfo.SeriesNumber=dicominfo.SeriesNumber+1000+SeriesNum;
dicominfo.SeriesInstanceUID=[dicominfo.SeriesInstanceUID(1:end-1), num2str(SeriesNum+1)];



mat2dicom_DCE(var,[imdir,filesep,Name,filesep],true,dicominfo);
disp([Name,' is saved as ',imdir, filesep ,Name])
SeriesNum=SeriesNum+1;
end
