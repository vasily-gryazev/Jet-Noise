%---------------------------------------------------%
%   Compute empirical temperature variance source 
%   for heated jet (Khavaran)
%---------------------------------------------------%
function ST = source_ST(js)

Dj = js.Dj; % jet diameter
NPR = js.NPR; % Nozzle Pressure Ratio
NTR = js.NTR; % Nozzle Temperature Ratio
LcD = js.LcD; % length of the potential core
X = js.X;
T = js.T;

Lc = LcD*Dj;
alpha = 0.2;
beta0 = 0.7;
delta = 1 + 1/(3*NPR);
C = (1-1/NTR)^delta/6;

XLc = min(X/Lc, 1);
beta = beta0 + (1-beta0)*XLc;

Tinf = min(T(end,:));
[T Tx Tr] = refine_data(js, T, Tinf);

ST = C*(abs(Tr)*Dj/Tinf).^alpha.*beta;
% ST = imgaussfilt(ST, 3);

end
