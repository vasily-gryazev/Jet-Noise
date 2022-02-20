warning('off','all') 

ks.par.phi = phi;
ks.data.Xg = ks.data.X(1,:)'; 
ks.data.Yg = ks.data.Y(:,1);
ks.par.modeN = 20;
ks.data.C2 = ks.data.C.^2;
ks.data.Source1 = ks.data.Kappa.^0.5/ks.data.Uj;
ks.data.Source2 = source2(ks.data);

t = linspace(0,1,nStg)';
Stg = 0.02*250.^t;
fs = {};
PSDg = [];
for i = 1:length(Stg)
    St = Stg(i);
    PSDg(i) = compute_PSD(ks, St);
    disp(i);
end

function [SPL PSD] = compute_PSD(ks, St)
data = ks.data;
x0 = max(data.Xg(1), 0.5*data.Dj);
Xg = linspace(x0, data.Xg(end), 10);

PSDg = [];
for i = 1:length(Xg)
    x = Xg(i);
    [P Pg] = compute_kav(ks, x, St);
    PSDg(i) = P;
end

PSD = trapz(Xg, PSDg);
const = constants();
a = data.Uj/data.Dj;
SPL = 10*log10(a*PSD/4e-10);
end

function [PSD Pg] = compute_kav(ks, x, St)
data = ks.data;
par = ks.par;

par = param_angle(par, par.phi);
par = param_St(data, par, St);

Rg = linspace(data.Yg(1), data.Yg(end), 200);
par = r_profile(data, par, x, Rg);

par.const.Cell = par.Cell;
par.const.Ctau = par.Ctau;

[PSD Pg] = r_SPD(par);
end

function par = r_profile(data, par, x0, Rg)
if nargin == 3, Rg = data.Yg; end
Rg = Rg(:);

rp = [];
rp.x0 = x0;
rp.Rg = Rg;

C = data_interp(data, 'C', x0, Rg);
rp.Cinf = C(end);
rp.k = par.omega/rp.Cinf;

rp.Ug = data_interp(data, 'U', x0, Rg);
Mg = rp.Ug/rp.Cinf;
Wg = 1 - Mg*par.cosphi;
rp.W2g = Wg.^2;
rp.W3g = Wg.^3;

Rhog = data_interp(data, 'Rho', x0, Rg);
rp.Rhoinf = Rhog(end);
rp.Phi2g = (Rhog/rp.Rhoinf).*rp.W2g;

rp.Mcg = data_interp(data, 'Uc', x0, Rg)/rp.Cinf;
rp.Kappag = data_interp(data, 'Kappa', x0, Rg);
rp.Epsg = data_interp(data, 'Eps', x0, Rg);
rp.Hrelg = data_interp(data, par.source, x0, Rg);
rp.Zetag = data_interp(data, 'Zeta', x0, Rg);
par.rp = rp;
end

function [iSPD SPD]= r_SPD(par)
rp = par.rp;
W2g = rp.W2g;
Mcg = rp.Mcg;
Phi2g = rp.Phi2g;
Kappag = rp.Kappag;
Epsg = rp.Epsg;
Hrelg = rp.Hrelg;
Rg = rp.Rg;

Rhoinf = rp.Rhoinf;
Cinf = rp.Cinf;
k = rp.k;
cosphi = par.cosphi;
sinphi = par.sinphi;

Ig = r_I1111(par);
Sg = r_green(par);

Qg = sqrt(Phi2g - cosphi^2);
R2g = (abs(cosphi^2 + Qg*sinphi)).^2;

c = Rhoinf^2*Ig*k^4.*W2g.*R2g.*Sg;
c = c./((4*pi*par.Rloc)^2*(1-Mcg*cosphi).^2);

SPDA = R2g.*c;
SPDB = W2g*(15/16).*(Cinf^2./Kappag).*Hrelg.^2.*c;

A = par.Aconst;
B = par.Bconst;
SPD = A*SPDA + B*SPDB;
iSPD = trapz(Rg, Rg.*SPD);
end

function Ig = r_I1111(par)
rp = par.rp;
const = par.const;
omegasg = (1-rp.Mcg*par.cosphi)*par.omega;
baru2g = rp.Kappag/1.5;

if strcmp(par.scales , 'k-eps')
    ellg = const.Cell*rp.Kappag.^1.5./rp.Epsg;
    tau0g = const.Ctau*rp.Kappag./rp.Epsg;
else
    ellg = const.Cell*sqrt(rp.Kappag)./rp.Zetag;
    tau0g = const.Ctau*1./rp.Zetag;
end

Hg = tau0g./(1+(omegasg.*tau0g/2).^2);
chig = (const.Cell/const.Ctau)*(par.omega.*tau0g).*(rp.Kappag.^0.5./343);
Ng = 5./(8*chig.^5).*(3*atan(chig)-chig.*(5*chig.^2+3)./(1+chig.^2).^2);
Ig = 4*ellg.^3/(5*pi^2).*baru2g.^2.*Hg.*Ng;
end


function [Sg Fg] = r_green(par)
rp = par.rp;
Rg = rp.Rg;
Phi2g = rp.Phi2g;
W3g = rp.W3g;
k = rp.k;
cosphi = par.cosphi;
sinphi = par.sinphi;
Ng = 0:(par.modeN-1);
NN = 20000;

Rinf = Rg(end); h = Rinf/NN;
rr = transpose(linspace(h, Rinf-h, NN-1));

[V V1] = spline_interp(Rg, Phi2g, rr); 

A = 1./rr - V1./V;
B0 = k^2*(V - cosphi^2);

aa = [1-A*h/2; 0; 0];
bb = [0; 0; 1+A*h/2];
cc0 = [0; B0*h^2-2; 0];
cc1 = [0; h^2./rr.^2; 0];
dd = zeros(NN+1, 1);
rr = [0; rr; Rinf];

kinf = k*sinphi;
Sg = zeros(length(Rg), 1);
if nargout == 2, Fg = zeros(length(Rg), length(Ng)); end

for i = 1:length(Ng)
    n = Ng(i);
    if n == 0
        eps = 1;
        alpha0 = 0; beta0 = 1;
        cc = cc0;
    else
        eps = 2*(-1i)^n;
        alpha0 = 1; beta0 = 0;
        cc = cc0 - n^2*cc1;
    end
    [betaR alphaR] = bessel_h(n, Rinf, kinf);
    gammaR = 2i*eps/(pi*Rinf);

    bb(2) = beta0;
    cc(1) = alpha0*h - beta0;
    aa(NN) = betaR;
    cc(NN+1) = alphaR*h - betaR;
    dd(NN+1) = gammaR*h;

    S = spdiags([aa cc bb], -1:1, NN+1, NN+1);
    gg = S\dd;

    ff = interp1(rr, gg, Rg)./W3g;
    if nargout == 2, Fg(:,i) = ff; end

    if i == 1, del = 2; else del = 1; end
    Sg(:) = Sg(:) + del*ff.*conj(ff);
end
end

function const = constants()
const = [];
const.gamma = 7/5;
const.R = 287;
const.Cv = const.R/(const.gamma-1);
const.Cp = const.gamma*const.Cv;
const.Pref2 = (2e-5)^2;
end

function [y y1] = bessel_j1(n, x, a)
x = a*x(:); y = besselj(n,x);
y1 = a/2*(besselj(n-1,x)-besselj(n+1,x));
end

function [y y1] = bessel_h(n, x, a)
x = a*x(:); y = besselh(n,x);
y1 = a/2*(besselh(n-1,x)-besselh(n+1,x));
end

function [sy sy1] = spline_interp(x, y, sx)
s = spline(x, y); sy = ppval(s, sx); 
s.coefs = s.coefs*diag(3:-1:1,1);
sy1 = ppval(s, sx); 
end

function V = data_interp(data, field, X, Y)
if nargin == 3, Y = data.Yg; end

if strcmp(field, 'C') 
    V = data_interp(data, 'C2', X, Y).^0.5;
elseif strcmp(field, 'Uc')
    a = data.aConv; b = data.bConv;
    V = a*data_interp(data, 'U', X, Y) + b*data.Uj;
elseif strcmp(field, 'Uc0') 
    V = 0*data_interp(data, 'U', X, Y);
elseif strcmp(field, 'U/Uj')
    V = data_interp(data, 'U', X, Y)/data.Uj;
else
    val = getfield(data, field);
    V = interp2(data.X, data.Y, val, X, Y);
end
end

function par = param_St(data, par, St)
if nargin == 2, St = 1; end

freq = (data.Uj/data.Dj)*St;
par = param_freq(data, par, freq);
end

function par = param_freq(data, par, freq)
if nargin == 1, freq = 1000; end

par.freq = freq; 
par.omega = 2*pi*freq; 
par.St = (data.Dj/data.Uj)*freq; 
end

function par = param_angle(par, phi)
if nargin == 1, phi = 90; end

par.phi = phi; 
phi = degtorad(phi);
par.cosphi = cos(phi);
par.sinphi = sin(phi);
end


