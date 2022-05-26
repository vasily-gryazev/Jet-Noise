%% Implementaion of Khavarn HOT JET NOISE model    
%   
% Input file should contain the variables in matrix form:
%      X     - axial coordinate 
%      Y     - radial coordinate
%      U     - axial velocity
%      Rho   - density 
%      T     - temperature
%      C     - speed of sound
%      Kappa - rurbulent kinetic energy
%      Eps   - turbulent dissipation rate
%      Zeta  - meanflow vorticity

%   Details can be found in:
%   [1] Khavaran, A., Kenzakowski, D.C., Mielke-Fagan, A.F., 
%   Hot jets and sources of jet noise, International Journal of Aeroacoustics, 
%   Vol. 9 pp. 491-532. https://doi.org/10.1260/1475-472X.9.4-5.491.  
%   
%   [2] Gryazev, V., Markesteijn, A.P., Karabasov, S.A.
%   On the robustness of reduced-order jet noise models, 
%   AIAA Journal (submitted), 2022.   
% 
%   Vasily Gryazev (v.gryazev@qmul.ac.uk) 
%   Last revision: May-2022

clear, clc, close all

%% Load meanflow data in matrix format and inspect dataset
%   The jet data stored in matrix form (200x500)
%   200 points in the radial and by 500 points in axial direction
%    
%   Please use 'permute' function to swap indices if necessary .
addpath('INPUT')

ks.data = load('SP7');
size(ks.data.X)

%% Define some parameters for SP7 jet 
%  Choose the definition of acoustic scales:
% 'k-eps' - turbulen kinetic energy + dissipation rate defined by eq.(6)
% 'k-zeta' - turbulen kinetic energy + dissipation rate defined by eq.(7)
%  and the corresponding calibarion constants ps.Cell and ps.Ctau 
%  in Appendix A
%  in AIAA paper: "On the robustness of reduced-order jet noise models"   

ks.data.name = 'SP7';           % jet name
ks.data.Uj   = 302;             % jet exit velocity                  
ks.data.Dj   = 0.0508;          % jet diameter
ks.data.NPR  = 2.861;           % Nozzle Pressure Ratio
ks.data.NTR  = 0.835;           % Nozzle Temperature Ratio
ks.data.LcD  = 5;               % length of the potential core in jet diameters (Dj)
ks.par.Rloc  = 100*ks.data.Dj;  % Observer location
 
% Select source model:
ks.par.source = 'Source1';      % Cold source model 
% ks.par.source = 'Source2';    % Hot source model

% Selecet acoustic scales definition
% Calibration constants see Table 5 in AIAA paper in Appendix A
% in AIAA paper: "On the robustness of reduced-order jet noise models"    

ks.par.scales = 'k-eps';        % TKE + dissipation rate 
% ks.par.scales = 'k-vort';     % TKE + meanflow vorticity
ks.par.Cell   = 0.62;           % cl
ks.par.Ctau   = 0.94;           % ctau
ks.par.Aconst = 9;              % Amplitude coefficient momentum source
ks.par.Bconst = 1;              % Amplitude coefficient enthalpy source
ks.data.aConv = 0;              % Convection velocity a
ks.data.bConv = 0.82;           % Convection velocity 

%% Finds acoustic integral and plot the spectrum
addpath('POST_PROC')

phi  = 90;   % Define observer angle
nStg = 5;    % Define number of frequencies (Strouhal number = f*Dj/Uj) 

run post_proc

resfile = sprintf('OUTPUT/%s_PSD_%d.mat',ks.data.name, phi);
save(resfile, '-struct', 'fs');

semilogx(fs.Stg, fs.PSDg,'r-o', 'linewidth', 2); 
xlabel('St'); ylabel('PSD')

%% Plot result vs experiment
%  Load experimental data and plot vs predictions 

ex = load(sprintf('EXP_DATA/EXP_%s_PSD_%d',ks.data.name, phi));
fs = load(sprintf('OUTPUT/%s_PSD_%d', ks.data.name, phi));

semilogx(ex.St, ex.PSD); hold on
semilogx(fs.Stg, fs.PSDg);
xlabel('St'); ylabel('PSD')










