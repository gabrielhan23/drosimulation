addpath(genpath('.'))
dir = "./results/#12 - comparison calculated";
mkdir(dir)

DRO = SimulatedDRO();
% #13 - pinn paper 
% DRO.flow = [0.01, 3.5];
% DRO.ps = [0.01, 4];
% DRO.vp = [0.01, .24];
% DRO.ve = [0.01, 0.6];

% #12 - comparison calculated
DRO.flow = [0.15, .5];
DRO.ps = [0.01, 0.1];
DRO.vp = [0.15, .3];
DRO.ve = [0.15, 0.5];

% #11 - best overall
% DRO.flow = [0.05, 1];
% DRO.ps = [0.01, 0.2];
% DRO.vp = [0.01, .4];
% DRO.ve = [0.01, 0.5];
simulate_make_dro(DRO, dir, (0:15:120), 15);
