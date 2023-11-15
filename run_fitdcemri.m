clc; close all; clear all;
% Load data with know ktrans and ve
% toi with ktrans=0.25 and ve=0.5
% rr with ktrans=0.1 and ve=0.1

load run_fitdcemri.mat;

%   1) Linear Reference Region Model
          pars1 = fitdcemri(toi,rr,time,'lsq');

%  2) Non-Linear Reference Region Model
x0=randn(3,1);
lb=zeros(3,1);
ub=10*ones(3,1);
          pars2 = fitdcemri(toi,rr,time,x0,lb,ub,'RRM');

%   3) Linear Tofts Model
          pars3 = fitdcemri(toi,Cp,time);

%   4) Non-linear Tofts Model
          pars4 = fitdcemri(toi,Cp,time,x0,lb,ub,'Tofts');
 
% Display RR parameters
          [pars1,pars2]
          
% Display Tofts parameters
          [pars3,pars4]