function gs = green_cr(js, gs)

Cinf = js.Cinf;
Dj = js.Dj;
inRg = js.Y(:,1);

gs.NN = 2000;
gs.ksmooth = ceil(gs.NN/50); 
gs.ksmooth1 = {'lowess', ceil(gs.NN/50)}; 
gs.ksmooth2 = {'lowess', ceil(gs.NN/100)}; 
gs.ksmoothG2 = ceil(gs.NN/10); 

gs.GF = {};

for iSt = 1:length(gs.Stg)
    St = gs.Stg(iSt);
    fprintf('St(%d:%d)=%1.2f\n', iSt, length(gs.Stg), St);

    freq = St*js.Uj/Dj;
    omega = 2*pi*freq;
    k = omega/js.Cinf;
    k1 = k*cosd(gs.phi);
    ks = k*sind(gs.phi);

    modemax = get_modemax(ks*Dj);

    Green(1:4,1:modemax+1) = {zeros(length(gs.Rg),length(gs.Xg))};

    for ix = 1:length(gs.Xg)
        xs = gs.Xg(ix);

        inUg = interp2(js.X, js.Y, js.U, xs, inRg);
        inCg = interp2(js.X, js.Y, js.C, xs, inRg);

        [Rg Ug C2g] = expand_grid(inRg, inUg, inCg, js.Cinf);

        NN = gs.NN;
        Rmax = Rg(end);
        RRg = linspace(0,Rmax,NN)';
        hh = RRg(2) - RRg(1);

        UUg = spline(Rg, Ug, RRg);
        CC2g = spline(Rg, C2g, RRg);

        [UUg UUrg UUrrg CC2g CC2rg] = diffusion(gs, RRg, UUg, CC2g, Dj);

        HHg = ones(NN,1);
        MMg = UUg/js.Cinf*cosd(gs.phi); 
        Mc = 0.99; 
        Icg = find(diff(MMg>Mc));
        if ~isempty(Icg) 
            delta = 0.05*Dj;
            idelta = ceil(delta/Rmax*NN);
            for ic = Icg(:)'
                i1 = max(ic-idelta,1);
                i2 = min(ic+idelta,NN);
                HHg(i1:i2) = HHg(i1:i2)*(1-1i*delta);
            end
            for i = 1:2 
                HHg = movmean(HHg, idelta);
            end
        end
        KK1g = HHg*k1;
        OOmega1g = omega - KK1g.*UUg;

        RR1g = [0; 1./RRg(2:end)];
        AAg = RR1g + CC2rg./CC2g;
        BBg = AAg.*UUrg + 2*KK1g./OOmega1g.*UUrg.^2 + UUrrg;
        BBg = OOmega1g.^2./CC2g - KK1g.^2 - KK1g./OOmega1g.*BBg;
        FFAg = RR1g - AAg;
        FFBg = ks^2 - BBg;

        Ampg = 1i*omega/(2*js.Rloc*js.Cinf^2)*exp(1i*(k*js.Rloc-KK1g*xs));

        for n = 0:modemax
            A = AAg(2:NN-1);
            B = BBg(2:NN-1);
            R1 = RR1g(2:NN-1);
            hN1 = besselh(n,ks*(Rmax-hh));
            hN = besselh(n,ks*Rmax);
            aa = [1-A*hh/2; hN; 0];
            bb = [0; 0; 1+A*hh/2];
            cc = [0; -2+(B-n^2*R1.^2)*hh^2; -hN1];

            Ang = Ampg*(-1i)^n;
            if n == 0
                cc(1) = 1;
                bb(2) = -1;
            else
                Ang = Ang*2;
                cc(1) = 1;
                bb(2) = 0;
            end

            PPig = besselj(n,ks*RRg);
            PPirg = ks*dbesselj(n,ks*RRg);
            FFg = FFAg.*PPirg + FFBg.*PPig;
            F = FFg(2:NN-1);
            ff = [0; F*hh^2; 0];

            SS = spdiags([aa cc bb], -1:1, NN, NN);
            PPg = SS\ff;

            PPg = PPg + PPig;
            PPrg = gradient(PPg, hh);

            PPg = PPg.*Ang;
            PPrg = PPrg.*Ang;

            G = {};
            G{1} = KK1g.*CC2g./OOmega1g.*PPg;
            G{2} = 1i./OOmega1g.*(CC2g.*PPrg-G{1}.*UUrg);
            if isfield(gs, 'ksmoothG2')
                G{2} = smooth(G{2}, gs.ksmoothG2, 'lowess');
            end
            PR1 = [PPrg(1); PPg(2:NN)./RRg(2:NN)];
            G{3} = -1i*n./OOmega1g.*CC2g.*PR1;
            G{4} = PPg;

            for m = 1:4
                Green{m,n+1}(:,ix) = interp1(RRg, G{m}, gs.Rg);
            end
        end
    end
    gs.GF{iSt} = Green;
end 

end

function modemax = get_modemax(ksDj)
    t = linspace(0,1,200);
    N = 15;
    modemax = N+1;
    for n = N:-1:0
        if max(abs(2*besselj(n,ksDj*t))) < 0.1
            modemax = modemax - 1;
        else
            break;
        end
    end
    modemax = min(modemax, 8);
end

function [Rg Ug C2g] = expand_grid(inRg, inUg, inCg, Cinf)
    N = length(inRg);
    M = ceil(0.35*N);
    Rg = zeros(N+M,1);
    Ug = zeros(N+M,1);
    C2g = zeros(N+M,1);

    Rg(1:N) = inRg;
    Ug(1:N) = inUg;
    C2g(1:N) = inCg.^2;

    h = Rg(N) - Rg(N-1);
    C2inf = Cinf^2;

    for i = N+1:length(Rg)
        Rg(i) = Rg(i-1) + h;
        e = exp(-((i-N)/10)^2);
        Ug(i) = Ug(N)*e;
        C2g(i) = C2inf + (C2g(N)-C2inf)*e;
    end
end

function [UUg UUrg UUrrg CC2g CC2rg] = diffusion(gs, RRg, UUg, CC2g, Dj)
    hh = RRg(2) - RRg(1);
    rad = Dj/2;

    UUg = smooth(UUg, gs.ksmooth, 'lowess');
    CC2g = smooth(CC2g, gs.ksmooth, 'lowess');
    UUg = max(UUg, 0);

    UUrg = [0; (UUg(3:end)-UUg(1:end-2))/(2*hh); 0];
    CC2rg = [0; (CC2g(3:end)-CC2g(1:end-2))/(2*hh); 0];

    I = find(RRg/rad>4);
    E = exp(-(RRg(I)/rad-4).^2);
    UUrg(I) = UUrg(I).*E;
    CC2rg(I) = CC2rg(I).*E;

    UUrg = smooth(UUrg, gs.ksmooth1{2}, 'lowess');
    CC2rg = smooth(CC2rg, gs.ksmooth1{2}, 'lowess');

    UUrrg = [0; (UUrg(3:end)-UUrg(1:end-2))/(2*hh); 0];
    UUrrg(1) = UUrrg(2);

    if strcmp(gs.ksmooth2{1},'lowess')
        UUrrg = smooth(UUrrg, gs.ksmooth2{2}, 'lowess');
    else
        for i = 1:gs.ksmooth2{2}
            D = 0.49*(UUrrg(1:end-4)-2*UUrrg(3:end-2)+UUrrg(5:end));
            UUrrg(3:end-2) = UUrrg(3:end-2) + D;
        end
    end
end

function jx = dbesselj(n,x)
    jx = 1/2*(besselj(n-1,x)-besselj(n+1,x));
end

function hx = dbesselh(n,x)
    hx = 1/2*(besselh(n-1,x)-besselh(n+1,x));
end
