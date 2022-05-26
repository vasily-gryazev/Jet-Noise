function ST = source2(data)

NPR = data.NPR;
NTR = data.NTR;
LcD = data.LcD;

Lc = LcD*data.Dj;
alpha = 0.2;
beta0 = 0.7;
delta = 1 + 1/(3*NPR);
C = (1-1/NTR)^delta/6;

ST = zeros(length(data.Yg), length(data.Xg));
for j = 1:length(data.Xg)
    x = data.Xg(j);
    t = x/Lc;
    if t <= 1
        beta = beta0 + (1-beta0)*t;
    else
        beta = 1;
    end
    Tg = data_interp(data, 'T', x);
    Tinf = Tg(end);
    [T T1] = spline_interp(data.Yg, Tg, data.Yg);
    for i = 1:length(data.Yg)
        y = data.Yg(i);
        dT = T1(i);
        ST(i,j) = C*(abs(dT)*data.Dj/Tinf)^alpha*beta;
    end
end

end

function [sy sy1] = spline_interp(x, y, sx)
% spline interpolation of y(x) and dy(x)/dx

s = spline(x, y); sy = ppval(s, sx); % y(sx)
s.coefs = s.coefs*diag(3:-1:1,1);
sy1 = ppval(s, sx); % dy(sx)/dx

end

function V = data_interp(data, field, X, Y)
% interpolate data field

% default Y profile
if nargin == 3, Y = data.Yg; end

if strcmp(field, 'C') % sound c
    V = data_interp(data, 'C2', X, Y).^0.5;
elseif strcmp(field, 'Uc')
    a = data.aConv; b = data.bConv;
    V = a*data_interp(data, 'U', X, Y) + b*data.Uj;
elseif strcmp(field, 'Uc0') % zero
    V = 0*data_interp(data, 'U', X, Y);
elseif strcmp(field, 'U/Uj')
    V = data_interp(data, 'U', X, Y)/data.Uj;
else
    val = getfield(data, field);
    V = interp2(data.X, data.Y, val, X, Y);
end

end



