function fs = post_proc(js, gs, ps)

ps.Nr = 150;
X = gs.X; 
R = gs.R;

if ps.Nr > length(gs.Rg)
    gs.Rg = linspace(gs.Rg(1),gs.Rg(end),ps.Nr)';
    [X R] = meshgrid(gs.Xg, gs.Rg);
end

U = refine_data(js, js.U, 0, X, R);
C = refine_data(js, js.C, js.Cinf, X, R);

Kappa = interp2(js.X, js.Y, js.Kappa, X, R);
Eps = interp2(js.X, js.Y, js.Eps, X, R);
Zeta = interp2(js.X, js.Y, js.Zeta, X, R); 

Rho = js.Rho;
Rhoinf = min(Rho(end,:));
Rho = refine_data(js, Rho, Rhoinf, X, R);
qhatsc2 = ps.A^2*(2/3*Rho.*Kappa).^2;

if strcmp(ps.scales, 'k-eps')
    ells = ps.Cell.*Kappa.^(3/2)./Eps;
    taus = ps.Ctau.*Kappa./Eps;
else
    ells = ps.Cell.*sqrt(Kappa)./Zeta;
    taus = ps.Ctau.*1./Zeta;
end

PSDg = zeros(size(gs.Stg));

for iSt = 1:length(gs.Stg)
    St = gs.Stg(iSt);

    freq = St*js.Uj/js.Dj;
    omega = 2*pi*freq;
    k = omega/js.Cinf;
    k1 = k*cosd(gs.phi);
    omega1 = omega - k1*U;
    W = 4*pi*(pi/log(2))^(3/2); 
    W = W.*qhatsc2.*ells.^3./taus;
    W = W.*exp(-(omega*ells./U).^2/(4*log(2)));
    W = W./(1+omega^2*taus.^2.*(1-U*cosd(gs.phi)/js.Cinf).^2);

    GF = gs.GF{iSt};
    totalS = zeros(size(R));

    for n = 0:size(GF,2)-1
        G4 = GF{4,n+1};
        if ps.Nr > length(gs.Rg)
            G4 = interp2(gs.X, gs.R, G4, X, R);
        end
        S = abs(G4).^2;
        if n == 0
            totalS = totalS + 2*S;
        else
            totalS = totalS + S;
        end
    end 
    totalS = totalS / 1000;
    
%     if gs.phi == 90
%         totalS = omega^2/(64*pi^4*js.Cinf^4*js.Rloc^2);
%     end
 
    P = W.*totalS;
    Pint = trapz(gs.Xg,P.*R,2);
    Pint = trapz(gs.Rg,Pint);
    PSDg(iSt) = Pint;
end

fs = {};
fs.Stg = gs.Stg;
fs.PSDg = 10*log10(4*pi*PSDg/(4e-10*(js.Dj/js.Uj)));

end