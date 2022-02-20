function fs = post_proc(js, gs, ps)

ps.Nr = 150; 
ps.r0cut = 0.1*js.Dj;
ps.epscut = 0.4*js.Dj; 
X = gs.X; 
R = gs.R;

if ps.Nr > length(gs.Rg)
    gs.Rg = linspace(gs.Rg(1),gs.Rg(end),ps.Nr)';
    [X R] = meshgrid(gs.Xg, gs.Rg);
end

[U Ux Ur] = refine_data(js, js.U, 0, X, R);
[V Vx Vr] = refine_data(js, js.V, 0, X, R);
C = refine_data(js, js.C, js.Cinf, X, R);

Kappa = interp2(js.X, js.Y, js.Kappa, X, R);
Eps = interp2(js.X, js.Y, js.Eps, X, R);
Zeta = interp2(js.X, js.Y, js.Zeta, X, R);

Rhoinf = min(js.Rho(end,:));
Rho = refine_data(js, js.Rho, Rhoinf, X, R);

SourceST = source_ST(js);         

R1 = 1./R; R1(1,:) = 0;
div = Ux + Vr + V.*R1;

ps.Celld = ps.Cell;
ps.Ctaud = ps.Ctau;
ps.Cellp = ps.Cell;            
ps.Cellpd = ps.Celld;

Cell = interp_lip(ps.rlip, ps.Cell, gs.Rg, gs.Xg);
Cellp = interp_lip(ps.rlip, ps.Cellp, gs.Rg, gs.Xg);
Ctau = interp_lip(ps.rlip, ps.Ctau, gs.Rg, gs.Xg);

Celld = interp_lip(ps.rlip, ps.Celld, gs.Rg, gs.Xg);
Cellpd = interp_lip(ps.rlip, ps.Cellpd, gs.Rg, gs.Xg);
Ctaud = interp_lip(ps.rlip, ps.Ctaud, gs.Rg, gs.Xg);

if strcmp(ps.scales, 'k-eps')
    A = (2/3*Rho.*Kappa).^2;
    ells = Cell.*Kappa.^(3/2)./Eps;
    ellsp = Cellp.*Kappa.^(3/2)./Eps;
    taus = Ctau.*Kappa./Eps;
    
    H = interp2(js.X, js.Y, SourceST, X, R);
    B = Rhoinf*js.Cinf^4*Rho.*Kappa.*H.^2;
    ellsd = Celld.*Kappa.^(3/2)./Eps;
    ellspd = Cellpd.*Kappa.^(3/2)./Eps;
    tausd = Ctaud.*Kappa./Eps;
else
    A = (2/3*Rho.*Kappa).^2;
    ells = Cell.*sqrt(Kappa)./Zeta;
    ellsp = Cellp.*sqrt(Kappa)./Zeta;
    taus = Ctau.*1./Zeta;

    H = interp2(js.X, js.Y, SourceST, X, R);
    B = Rhoinf*js.Cinf^4*Rho.*Kappa.*H.^2;
    ellsd = Celld.*sqrt(Kappa)./Zeta;
    ellspd = Cellpd.*sqrt(Kappa)./Zeta;
    tausd = Ctaud.*1./Zeta;
end

C1111 = interp_lip(ps.rlip, ps.C1111, gs.Rg, gs.Xg);
C2222 = interp_lip(ps.rlip, ps.C2222, gs.Rg, gs.Xg);
C3333 = interp_lip(ps.rlip, ps.C3333, gs.Rg, gs.Xg);
C1212 = interp_lip(ps.rlip, ps.C1212, gs.Rg, gs.Xg);
C2323 = interp_lip(ps.rlip, ps.C2323, gs.Rg, gs.Xg);
C1313 = interp_lip(ps.rlip, ps.C1313, gs.Rg, gs.Xg);
C1122 = interp_lip(ps.rlip, ps.C1122, gs.Rg, gs.Xg);
C1133 = interp_lip(ps.rlip, ps.C1133, gs.Rg, gs.Xg);
D11 = interp_lip(ps.rlip, ps.D11, gs.Rg, gs.Xg);
D22 = interp_lip(ps.rlip, ps.D22, gs.Rg, gs.Xg);
D33 = interp_lip(ps.rlip, ps.D33, gs.Rg, gs.Xg);
D12 = interp_lip(ps.rlip, ps.D12, gs.Rg, gs.Xg);

PSD = {}; PSD(1:length(gs.Stg),1:2) = {zeros(length(gs.Rg),1)};
PSDx = {}; PSDx(1:length(gs.Stg),1:2) = {zeros(length(gs.Xg),1)};

for iSt = 1:length(gs.Stg)
    St = gs.Stg(iSt);

    freq = St*js.Uj/js.Dj;
    omega = 2*pi*freq;
    k = omega/js.Cinf;
    k1 = k*cosd(gs.phi);
    omega1 = omega - k1*U;
    
    W = (pi/log(2))^(3/2);
    W = W.*ells.*ellsp.^2.*taus;
    W = W./(1+(omega1.*taus).^2);
    W = W.*exp(-(omega*ells./U).^2/(4*log(2)));

    Wd = 2*pi*(pi/log(2))^(3/2);
    Wd = Wd.*ellsd.*ellspd.^2.*tausd;
    Wd = Wd./(1+(omega1.*tausd).^2);
    Wd = Wd.*exp(-(omega*ellsd./U).^2/(4*log(2)));
    
    GF = gs.GF{iSt};
    totalSA = {}; totalSA(1:2) = {zeros(size(R))}; % cold
    totalSB = {}; totalSB(1:2) = {zeros(size(R))}; % hot

    for n = 0:size(GF,2)-1
        G = {}; Gx = {}; Gr = {};
        for m = 1:4
            G{m} = GF{m,n+1};
            [Gx{m} Gr{m}] = gradient(G{m}, gs.Xg(2)-gs.Xg(1), gs.Rg(2)-gs.Rg(1));
            if ps.Nr > length(gs.Rg)
                G{m} = interp2(gs.X, gs.R, G{m}, X, R);
                Gx{m} = interp2(gs.X, gs.R, Gx{m}, X, R);
                Gr{m} = interp2(gs.X, gs.R, Gr{m}, X, R);
            end
        end

        ga = 1.4;
        G{4} = G{4}*(ga-1);
        Gx{4} = Gx{4}*(ga-1);
        Gr{4} = Gr{4}*(ga-1);

        Gammag = U.*Gx{4} + V.*Gr{4} - (ga-1)*div.*G{4} + 1i*omega*G{4};

        Ia11 = Gx{1} - (Ux.*G{4}+U.*Gx{4}) + Gammag/2;
        Ia12 = Gx{2} - (Vx.*G{4}+V.*Gx{4});
        Ib13 = Gx{3};
        Ia21 = Gr{1} - (Ur.*G{4}+U.*Gr{4});
        Vr_bett = -Gammag/(ga-1) - Gx{1} - R1.*(n*G{3}+G{2});
        Ia22 = Vr_bett - (Vr.*G{4}+V.*Gr{4}) + Gammag/2;
        Ib23 = Gr{3};
        Ib31 = -n*R1.*(G{1}-U.*G{4});
        Ib32 = -R1.*(n*G{2}+G{3}-n*V.*G{4});
        Ia33 = R1.*(n*G{3}+G{2}-V.*G{4}) + Gammag/2;

        Ja1 = Gx{4};
        Ja2 = Gr{4};
        Jb3 = -n*R1.*G{4};

        SA = ...
              C1111.*real(Ia11.*conj(Ia11)) + ...
            2*C1122.*real(Ia11.*conj(Ia22)) + ...
            2*C1133.*real(Ia11.*conj(Ia33)) + ...
              C1212.*real(Ia12.*conj(Ia12)) + ...
              C1212.*real(Ia21.*conj(Ia21)) + ...
            2*C1212.*real(Ia21.*conj(Ia12)) + ...
              C2323.*real(Ib23.*conj(Ib23)) + ...
              C1313.*real(Ib13.*conj(Ib13)) + ...
              C1313.*real(Ib31.*conj(Ib31)) + ...
            2*C1313.*real(Ib31.*conj(Ib13)) + ...
              C2222.*real(Ia22.*conj(Ia22)) + ...
              C3333.*real(Ia33.*conj(Ia33));

        SB = ...
              D11.*real(Ja1.*conj(Ja1)) + ...
              D22.*real(Ja2.*conj(Ja2)) + ...
              D33.*real(Jb3.*conj(Jb3)) + ...
            2*D12.*real(Ja1.*conj(Ja2));

        if n == 0
            totalSA{1} = totalSA{1} + 2*SA;
            totalSB{1} = totalSB{1} + 2*SB;
        else
            totalSA{2} = totalSA{2} + SA;
            totalSB{2} = totalSB{2} + SB;
        end
    end 

    for m = 1:2
        PA = W.*A.*totalSA{m}/(2*pi);
        PB = Wd.*B.*totalSB{m}/(2*pi);
        
        P = ps.Acold*PA + ps.Bhot*PB;
        Pint = trapz(gs.Xg,P.*R,2);
        PSD{iSt,m} = Pint(:);
        
        Pintx = trapz(gs.Rg,P.*R);
        PSDx{iSt,m} = Pintx(:);
    end
end 

PSDg = zeros(length(gs.Stg),1);
r0 = ps.r0cut; eps = ps.epscut;
ra = max(r0-eps/2,gs.Rg(1));
rb = min(r0+eps/2,gs.Rg(end));
NR = length(gs.Rg)/gs.Rg(end);
ia = max(ceil(ra*NR),1);
ib = min(ceil(rb*NR),length(gs.Rg));
ii = ia:ib;

for iSt = 1:length(gs.Stg)
    if ib > ia
        y = PSD{iSt,1};
        y(ii) = y(ia) + (y(ib)-y(ia))/(ib-ia)*(ii-ia);
        PSD{iSt,1} = y;
    end
    P = PSD{iSt,1} + PSD{iSt,2};
    PSDg(iSt) = trapz(gs.Rg,P);
end

a = js.Uj/js.Dj/4e-10;
PdBg = 10*log10(a*PSDg);

fs = {};                                                                  
fs.phi = gs.phi; 
fs.Stg = gs.Stg; 
fs.PSDg = PdBg;

matfile = sprintf('OUTPUT/%s_PSD_%d.mat',js.name, gs.phi);
save(matfile, '-struct', 'fs');
end

