% Leib's critical factor
L1111 = ones(size(X));
return

Nx = size(omega1,2);
Rcg = zeros(Nx,1);
RRg = linspace(Rg(1),Rg(end),100*length(Rg));
for ix = 1:Nx
    OOg = spline(Rg, omega1(:,ix), RRg);
    IIc = find(diff(OOg>0));
    if ~isempty(IIc)
        Rcg(ix) = RRg(IIc(end));
    end
end

if max(Rcg) == 0
    return
end

Rcg = movmean(Rcg, 5);
Rcxg = gradient(Rcg, Xg(2)-Xg(1));

eps = 0.5*Dj;
for ix = 1:Nx
    rc = Rcg(ix);
    if rc > 0
        rcx = Rcxg(ix);
        Uc = interp1(Rg, U(:,ix), rc);
        Vc = interp1(Rg, V(:,ix), rc);
        Ucr = interp1(Rg, Ur(:,ix), rc);
        alpha2 = 1i*k1*Ucr/(2*(Vc-Uc*rcx/eps));
        y2 = (Rg-rc).^2;
        z2 = alpha2.^y2;
        % L = abs(z2).*abs(igamma(1/2,z2)).^2;
        L = abs(z2).*abs(gammainc(1/2,abs(z2))).^2;
        L1111(:,ix) = L;
    end
end

if 0
    figure
    surf(X, R, L1111, 'edgecolor', 'none');
    colorbar; view(0,90);
    xlabel('X/Dj'); ylabel('R/Dj');
    axis tight
end
