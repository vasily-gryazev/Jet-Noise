% KHAVARAN HOT JET NOISE SOLVER    
%   
%   Details can be found in:
%   [1] Khavaran, A., Kenzakowski, D.C., Mielke-Fagan, A.F., 
%   Hot jets and sources of jet noise, International Journal of Aeroacoustics, 
%   Vol. 9 pp. 491-532. https://doi.org/10.1260/1475-472X.9.4-5.491.  
%   
%   [2] Gryazev, V., Markesteijn, A.P., Karabasov, S.A.
%   On the robustness of reduced-order jet noise models, 
%   AIAA Journal (submitted), 2022.   

clear, clc, close all

%% Load jet data
addpath('INPUT')
jet = 'SP7';
ks.data = load(jet);

ks.data.Uj = 302;
ks.data.Dj = 0.0508;
ks.data.NPR = 2.861;
ks.data.NTR = 0.835;
ks.data.LcD = 6;
ks.par.Rloc = 100*ks.data.Dj;

% Select source model:
ks.par.source = 'Source1';      % cold source  
% ks.par.source = 'Source2';    % hot source

% Selecet scales definition
ks.par.scales = 'k-eps';        % TKE + dissipation rate
% ks.par.scales = 'k-vort';     % TKE + meanflow vorticity

ks.par.Cell = 0.62;
ks.par.Ctau = 0.94;
ks.par.Aconst = 9;
ks.par.Bconst = 1;
ks.data.aConv = 0;
ks.data.bConv = 0.82;

%% Run 
addpath('POST_PROC')

phi = 30;   % observer angle
nStg = 10;  % generate n St numbers

run post_proc

fs = {};                                                                  
fs.phi = phi; fs.Stg = Stg; fs.PSDg = PSDg;

resfile = sprintf('OUTPUT/%s_PSD_%d.mat',jet, phi);
save(resfile, '-struct', 'fs');

%% Plot result vs experiment
jet = 'SP7';
phi = 30;

ex = load(sprintf('EXP_DATA/EXP_%s_PSD_%d',jet, phi));
fs = load(sprintf('OUTPUT/%s_PSD_%d', jet, phi));

semilogx(ex.St, ex.PSD); hold on
semilogx(fs.Stg, fs.PSDg);
xlabel('St'); ylabel('PSD')










