%% Impementaion of TAM and Auriault model 
%   
% Instructions to Matlab code.
% 
% The code has several iterations (steps) to calculate the spectrum:
% 
% Everything can be done via main routine
% step by step 3 structures are created to calculate the spectrum:
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
%   [1] Tam, C. K. W., Auriault, L., 
%   Jet Mixing Noise from Fine-Scale Turbulence, 
%   AIAA Journal, Vol. 37, No. 2, pp.145-153, 1999.  
%   
%   [2] Gryazev, V., Markesteijn, A.P., Karabasov, S.A.
%   On the robustness of reduced-order jet noise models, 
%   AIAA Journal (submitted), 2022.   

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

%% Calculate Green's function 
%  At the first run calculate Green's function 
%  which is then store in GREEN folder in accordance with
%  jet name and observer angle 
%  for example GREEN/SP3_green_90.
% 
%  In the next runs Green's function .mat file can be called from GREEN folder,
%  This is done to save time if model recalibration is required 
%  which is explained in the next section  
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

%%% choose the definition of th scales and the corresponding constants
ps.scales = 'k-eps'; ps.Cell = 0.265; ps.Ctau = 0.220; ps.A = 0.755;       
% ps.scales = 'k-vort';  ps.Cell = 0.950; ps.Ctau = 0.785; ps.A = 0.755;            

fs = post_proc(js, gs, ps); 

semilogx(fs.Stg, fs.PSDg,'r-o', 'linewidth', 2); 
xlabel('St'); ylabel('PSD')

resfile = sprintf('OUTPUT/%s_PSD_%d.mat',js.name, gs.phi);
save(resfile, '-struct', 'fs');

%% Plot result vs experiment
%  Load experimental data and plot vs predictions 

ex = load(sprintf('EXP_DATA/EXP_%s_PSD_%d',js.name, gs.phi));
fs = load(sprintf('OUTPUT/%s_PSD_%d', js.name, gs.phi));

semilogx(ex.St, ex.PSD, 'k:'); hold on
semilogx(fs.Stg, fs.PSDg,'r-o',  'linewidth', 1);
xlabel('St'); ylabel('PSD')


