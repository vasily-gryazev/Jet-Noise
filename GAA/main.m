%% GAA solver
%  
%   Details can be found in:
%   [1] Gryazev, V., Markesteijn, A.P., Karabasov, S.A.
%   Generalised acoustic analogy modelling of hot jet noise, 
%   AIAA Journal, Nov, 2021. https://doi.org/10.2514/1.J060896  
% 
%   [2] Gryazev, V., Markesteijn, A.P., Karabasov, S.A.
%   On the robustness of reduced-order jet noise models, 
%   AIAA Journal (submitted), 2022.   

clear, clc, close all

%% Load input meanflow data
addpath('INPUT')

jet  = 'SP3';
% jet  = 'SP7';

js = load(jet); 
js.name = jet;   

%% Green's function 
addpath('GREEN')

gs = {}; 
gs.phi = 30;     % observer angle
gs.nStg = 10;    % n of Strouhal number 

%%% calculate Green's funciton ~5-10 min
% gs = green_to_mat(js, gs);   
%%% or load existing Green's function
gs = load(sprintf('GREEN/%s_green_%d',js.name,gs.phi));         

%% Post-processing
addpath('POST_PROC')

ps = {};
ps.rlip = 0.5*js.Dj; 
ps.Acold = 0.25;           ps.Bhot = 0;
ps.C1111 = 1;              ps.D11 = 1;
ps.C2222 = 0.355;          ps.D22 = 0.55;
ps.C3333 = 0.360;          ps.D33 = 0.46;
ps.C1212 = 0.327;          ps.D12 = 0.01;
ps.C1313 = 0.326;
ps.C2323 = 0.180;
ps.C1122 = 0;
ps.C1133 = 0;

%%% choose the definition of th scales and the corresponding constants
% ps.scales = 'k-eps'; ps.Cell = 0.696; ps.Ctau = 0.551;    % TKE + dissipation rate           
ps.scales = 'k-zeta'; ps.Cell = 2.196; ps.Ctau = 1.333; % TKE + meanflow vorticity          

fs = post_proc(js, gs, ps); 
semilogx(fs.Stg, fs.PSDg,'r-o', 'linewidth', 2); 
xlabel('St'); ylabel('PSD')

resfile = sprintf('OUTPUT/%s_PSD_%d.mat',js.name, gs.phi);
save(resfile, '-struct', 'fs');

%%% Plot result vs experiment
ex = load(sprintf('EXP_DATA/EXP_%s_PSD_%d',jet, gs.phi));
fs = load(sprintf('OUTPUT/%s_PSD_%d', jet, gs.phi));

semilogx(ex.St, ex.PSD); hold on
semilogx(fs.Stg, fs.PSDg);
xlabel('St'); ylabel('PSD')











