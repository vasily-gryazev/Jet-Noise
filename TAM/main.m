%% TAM solver
%   
%   Details can be found in:
%   [1] Tam, C. K. W., Auriault, L., 
%   Jet Mixing Noise from Fine-Scale Turbulence, 
%   AIAA Journal, Vol. 37, No. 2, pp.145-153, 1999.  
%   
%   [2] Gryazev, V., Markesteijn, A.P., Karabasov, S.A.
%   On the robustness of reduced-order jet noise models, 
%   AIAA Journal (submitted), 2022.   

clear, clc, close all

%% Load data
addpath('INPUT')

jet  = 'SP7';  
js = load(jet); 
js.name = jet;   

%% Green's function
addpath('GREEN')

gs = {}; 
gs.phi = 90;
gs.nStg = 10; 

%%% calculate Green's funciton ~5-10 min
gs = green_to_mat(js, gs);   
%%% or load existing Green's function
% gs = load(sprintf('GREEN/%s_green_%d',js.name,gs.phi));     

%% 
addpath('POST_PROC')

ps = {};

%%% choose the definition of th scales and the corresponding constants
ps.scales = 'k-eps'; ps.Cell = 0.265; ps.Ctau = 0.220; ps.A = 0.755;        % TKE + dissipation rate   
% ps.scales = 'k-vort';  ps.Cell = 0.950; ps.Ctau = 0.785; ps.A = 0.755;    % TKE + meanflow vorticity          

fs = post_proc(js, gs, ps); 
semilogx(fs.Stg, fs.PSDg,'r-o', 'linewidth', 2); 
xlabel('St'); ylabel('PSD')

resfile = sprintf('OUTPUT/%s_PSD_%d.mat',js.name, gs.phi);
save(resfile, '-struct', 'fs');

%% Plot result vs experiment
jet = 'SP7';
phi = 90;

ex = load(sprintf('EXP_DATA/EXP_%s_PSD_%d',jet, phi));
fs = load(sprintf('OUTPUT/%s_PSD_%d', jet, phi));

semilogx(ex.St, ex.PSD); hold on
semilogx(fs.Stg, fs.PSDg);
xlabel('St'); ylabel('PSD')


