function Y = interp_lip(rlip, ylip, Rg, Xg)

Ri = [Rg(1); rlip(:); Rg(end)];
Yi = [ylip(1); ylip(:); ylip(end)];
Yg = interp1(Ri, Yi, Rg);
for i = 1:1
    Yg = movmean(Yg, 5);
end
Y = repmat(Yg, 1, length(Xg));

end
