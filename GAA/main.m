%% Implementaion of Generalised Acoustic Analogy solver
%  
% Instructions to Matlab code.
% 
% The code has several iterations (steps) to calculate the spectrum:
% 
% Everything can be done via main routine
% step by step  3 structures are created to calculate the spectrum:
% 
% 1. js = {} -- jet structure (Test datase example provided in the INPUT folder)
% 2. gs = {} -- green structure (Green's funcition)
% 3. ps = {} -- calibration parameters structure (cl, ctau, alpha_ijkl)
% 
% In some details:
% 1. Input is meanflowfields solution on the uniform equally spaced grid. 
% 	 Grid dimensions: 
% 	 500 points in x-direction (20Dj), 200 points in y-direction (5Dj).
% 
%    Input file should contain the variables in matrix form:
%      X     - axial coordinate 
%      Y     - radial coordinate
%      U     - axial velocity
%      V     - radial velocity
%      Rho   - density 
%      T     - temperature
%      C     - speed of sound
%      Kappa - rurbulent kinetic energy
%      Eps   - turbulent dissipation rate
%      Zeta  - meanflow vorticity
% 	
%   2. Calculating Greens Function and saving to .mat file:
% 	Calling Jet Structure (js) to compute
% 	solution for locally-parallel Green's function
% 	for observer angles phi = 90,60,30 and write
% 	the components of the vector adjoint
% 	Green's function (u, v, w, p) to Green Structure (gs)
% 	and saving it to .mat file.
% 		 
%   The resulting mat file is then stored in GREEN folder.
% 
%   3. Calculating the spectrum for each 
% 	Calling Jet Structure (js) and Green Stucture (gs)
% 	to calculate the final spectrum.
% 	ps. structure includes all calibration parameters.
% 
%   Once jet data is converted to mat file and the Greens function is found,
%   you can simply call these 2 matfiles from RANS and GREEN folders
%   This significantly reduce the time spent on model calibration. 
% 
%   Details can be found in:
%   [1] Gryazev, V., Markesteijn, A.P., Karabasov, S.A.
%   Generalised acoustic analogy modelling of hot jet noise, 
%   AIAA Journal, Nov, 2021. https://doi.org/10.2514/1.J060896  
% 
%   [2] Gryazev, V., Markesteijn, A.P., Karabasov, S.A.
%   On the robustness of reduced-order jet noise models, 
%   AIAA Journal (submitted), 2022.   


% Vasily Gryazev (v.gryazev@qmul.ac.uk)  
% Last revision: May-2022

clear, clc, close all

%% Load meanflow data in matrix format and inspect dataset
%   The jet data stored in matrix form (200x500)
%   200 points in the radial and by 500 points in axial direction
%    
%   Please use 'permute' function to swap indices if necessary .
addpath('INPUT'); 

js = load('SP7');   

size(js.X)

%% Define some parameters for SP7 jet 
js.name    = 'SP7';          % jet name
js.Dj      = 0.0508;         % jet diameter
js.Rloc    = 5.08;           % Observer location in meters
js.Cinf    = 340;            % Ambient Speed of sound (m/s) 

js.Uj      = 302;            % jet exit velocity
js.NPR     = 2.861;          % Nozzle Pressure Ratio (
js.NTR     = 0.835;          % Nozzle temperature ratio
js.LcD     = 5;              % length of the potential core in jet diametrs (Dj)

%% Calculate Green's function 
%  At the first run calculate Green's function 
%  which is then store in GREEN folder in accordance with
%  jet name and observer angle 
%  for example GREEN/SP3_green_90.
% 
%  In the next runs Green's function .mat file can be called from GREEN folder,
%  This is done to save time if model recalibration is required 
%  which is explained in the next section  
%
addpath('GREEN')

gs        = {}; 
gs.phi    = 90;              % Define observer angle
gs.nStg   = 5;               % Define number of frequencies (Strouhal number = f*Dj/Uj) 

%%%% calculate Green's funciton ~5-10 min
gs = green_to_mat(js, gs);   
%%%% load Green's function from GREEN folder  
% gs = load(sprintf('GREEN/%s_green_%d',js.name, gs.phi));     

%% Calculate Acoustic Spectrum
%  Finds acoustic integral and plot the spectrum
%  Choose the definition of acoustic scales:
% 'k-eps' - turbulen kinetic energy + dissipation rate defined by eq.(6)
% 'k-zeta' - turbulen kinetic energy + dissipation rate defined by eq.(7)
%  and the corresponding calibarion constants ps.Cell and ps.Ctau 
%  in Appendix A
%  in AIAA paper: "On the robustness of reduced-order jet noise models"    
%   
%  Result is saved to OUTPUT folder

addpath('POST_PROC')

ps = {};
ps.COLD_SOURCE = 1;        % Cold source (1 always)
ps.HOT_SOURCE  = 0;        % Hot source (change to 1 to activate)

ps.scales      = 'k-eps'; ps.Cell = 0.696; ps.Ctau = 0.551;             
% ps.scales      = 'k-zeta'; ps.Cell = 2.196; ps.Ctau = 1.333;    

ps.rlip        = 0.5*js.Dj; % lipline location 
ps.C1111       = 1;         % Nondimensional amplitude parametrs of GAA model            
ps.C2222       = 0.355;     % for fluctuating momentum source   
ps.C3333       = 0.360;     % see Table 4 in AIAA paper
ps.C1212       = 0.327;     
ps.C1313       = 0.326;      
ps.C2323       = 0.180;
ps.C1122       = 0;
ps.C1133       = 0;

ps.D11         = 1;         % Nondimensional amplitude parametrs of GAA model  
ps.D22         = 0.55;      % for the fluctuating enthalpy source term
ps.D33         = 0.46;      % see Table 5 in AIAA paper
ps.D12         = 0.01;      

fs = post_proc(js, gs, ps); 

semilogx(fs.Stg, fs.PSDg,'r-o', 'linewidth', 1);
xlabel('St');ylabel('PSD')

resfile = sprintf('OUTPUT/%s_PSD_%d.mat',js.name, gs.phi);
save(resfile, '-struct', 'fs');   

%% Plot result vs experiment
%  Load experimental data and plot vs predictions 

ex = load(sprintf('EXP_DATA/EXP_%s_PSD_%d',js.name, gs.phi));
fs = load(sprintf('OUTPUT/%s_PSD_%d', js.name, gs.phi));

semilogx(ex.St, ex.PSD, 'k:'); hold on
semilogx(fs.Stg, fs.PSDg,'r-o',  'linewidth', 1);
xlabel('St'); ylabel('PSD')













