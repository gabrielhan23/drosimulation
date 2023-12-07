%%Plot from excel sheet

for n=1:4;
    DCEdata(n,1,1:4)=fitdcemri(Data.R1.SofiD5(n,:)-0.76923,Data.R1.SofiD5(3,:)-0.76923,Data.time,'lsq');
    DCEdata(n,2,1:4)=fitdcemri(Data.R1.ParisD8(n,:)-0.76923,Data.R1.ParisD8(3,:)-0.76923,Data.time,'lsq');
    DCEdata(n,3,1:4)=fitdcemri(Data.R1.LisbonD6(n,:)-0.76923,Data.R1.LisbonD6(3,:)-0.76923,Data.time,'lsq');
end
%%
SubjectN=[1,2,3];
ParaN=[1,2,4];

figure;scatter(Ktrans(1,SubjectN),Ve(1,SubjectN),'b')%Remote
hold on;scatter(Ktrans(2,SubjectN),Ve(2,SubjectN),'r')%MVO
hold on;scatter(Ktrans(4,SubjectN),Ve(4,SubjectN),'g')%MI

%%
hold on;scatter(Ktrans(ParaN,1),Ve(ParaN,1),'v')%sofi
%hold on;scatter(Ktrans(ParaN,2),Ve(ParaN,2),'x')%Patis
hold on;scatter(Ktrans(ParaN,3),Ve(ParaN,3),'*')%Lisbon

%% Prediction
Ptime=0:1:1800;
%
rr=Data.R1.SofiD5(3,:)-0.76923;
PCp=fit(time',rr','exp2')
Cp=PCp(Ptime');
figure;plot(PCp,time',rr');
 [Ktrans,Kep,KepRR,Rsq]=fitdcemri(Data.R1.SofiD5(n,:)-0.76923,Data.R1.SofiD5(3,:)-0.76923,Data.time,'lsq');
Cp='a1*exp(-m1*x)+a2*exp(-m2*x)'
fitfun=fittype(@(a1,a2,m1,m2,x) a1*exp(-m1*x)+a2*exp(-m2*x))
x0=[5 5 0 0];
Pcp=fit(time',rr'/max(rr),fitfun,'StartPoint',x0)
%%LM model eq21
 Ct='DKtrans*(a1/(m1-kep)*(exp(-(Ktrans/ve)*x)-exp(-m1*x))+a2/(m2-kep)*(exp(-(Ktrans/ve)*x)-exp(-m2*x)))'
 %Ct=['DKtrans*(',num2str(Pcp.a1),'/(',num2str(Pcp.m1),'-kep)*(exp(-(Ktrans/ve)*x)-exp(-',num2str(Pcp.m1),'*x))+'...
 %    ,num2str(Pcp.a2),'/(',num2str(Pcp.m2),'-kep)*(exp(-(Ktrans/ve)*x)-exp(-',num2str(Pcp.m2),'*x)))'];
 %Ctt=fit(time',Data.R1.SofiD5(n,:)'-0.76923,'exp3')

 Ctt=fit(time',Data.R1.SofiD5(n,:)'-0.76923,Ct)
 figure;plot(Ctt,time,Data.R1.SofiD5(n,:)'-0.76923)

%%
x0=randn(3,1);
lb=zeros(3,1);
ub=10*ones(3,1);
for n=1:4;
    DCEdata(n,1,1:4)=fitdcemri(Data.R1.SofiD5(n,:)-0.76923,Data.R1.SofiD5(3,:)-0.76923,Data.time,'lsq');
    DCEdata(n,2,1:4)=fitdcemri(Data.R1.ParisD8(n,:)-0.76923,Data.R1.ParisD8(3,:)-0.76923,Data.time,'lsq');
    DCEdata(n,3,1:4)=fitdcemri(Data.R1.LisbonD6(n,:)-0.76923,Data.R1.LisbonD6(3,:)-0.76923,Data.time,'lsq');
end
%%
SubjectN=[1,2,3];
ParaN=[1,2,4];

figure;scatter(Ktrans(1,SubjectN),Ve(1,SubjectN),'b')%Remote
hold on;scatter(Ktrans(2,SubjectN),Ve(2,SubjectN),'r')%MVO
hold on;scatter(Ktrans(4,SubjectN),Ve(4,SubjectN),'g')%MI

%%
hold on;scatter(Ktrans(ParaN,1),Ve(ParaN,1),'v')%sofi
%hold on;scatter(Ktrans(ParaN,2),Ve(ParaN,2),'x')%Patis
hold on;scatter(Ktrans(ParaN,3),Ve(ParaN,3),'*')%Lisbon

function Ct =LMmodel(Cp,Ktrans,Ve,Ptime)
Ct=conv(Cp,Ktrans*exp(-Ktrans/Ve)*Ptime)
end