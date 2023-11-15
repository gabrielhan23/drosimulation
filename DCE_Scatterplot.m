%scatterplot
function []= DCE_Scatterplot(fitresultsDCEetk,RemoteMask,MIMask,MVOMask)
% Remote
temp3=fitresultsDCEetk(:,:,3);
temp3=temp3(RemoteMask==1);
temp2=fitresultsDCEetk(:,:,2);
temp2=temp2(RemoteMask==1);
temp1=fitresultsDCEetk(:,:,1);
temp1=temp1(RemoteMask==1);
figure;scatter3(temp1(:),temp2(:),temp3(:),'o','r')

% MI
temp3=fitresultsDCEetk(:,:,3);
temp3=temp3(MIMask==1);
temp2=fitresultsDCEetk(:,:,2);
temp2=temp2(MIMask==1);
temp1=fitresultsDCEetk(:,:,1);
temp1=temp1(MIMask==1);
hold on;scatter3(temp1(:),temp2(:),temp3(:),'x','b')


% MVO
temp3=fitresultsDCEetk(:,:,3);
temp3=temp3(MVOMask==1);
temp2=fitresultsDCEetk(:,:,2);
temp2=temp2(MVOMask==1);
temp1=fitresultsDCEetk(:,:,1);
temp1=temp1(MVOMask==1);
hold on;scatter3(temp1(:),temp2(:),temp3(:),'v','g')

scalemask=(RemoteMask==1)|(MIMask==1)|(MVOMask==1);
temp3=fitresultsDCEetk(:,:,3);

temp2=fitresultsDCEetk(:,:,2);

temp1=fitresultsDCEetk(:,:,1);
scalemask(isnan(temp1)|isnan(temp2)|isnan(temp3))=0;

Mn=[mean(temp1(scalemask)),mean(temp2(scalemask)),mean(temp3(scalemask))]
SD=[std(temp1(scalemask)),std(temp2(scalemask)),std(temp3(scalemask))]
xlim([Mn(1)-2*SD(1) Mn(1)+2*SD(1)])
ylim([Mn(2)-2*SD(2) Mn(2)+2*SD(2)])
zlim([Mn(3)-2*SD(3) Mn(3)+2*SD(3)])