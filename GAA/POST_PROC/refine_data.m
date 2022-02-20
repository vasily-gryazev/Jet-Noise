function [F Fx Fy] = refine(js, F, Finf, Xq, Yq)

X = js.X;
Y = js.Y;
Dj = js.Dj;

I = find(Y>2*Dj);
E = exp(-(Y(I)/Dj-2).^2);
F(I) = Finf + (F(I)-Finf).*E;

if nargout > 1
    [Fx Fy] = gradient(F, X(1,2)-X(1,1), Y(2,1)-Y(1,1));
    Fy(I) = Fy(I).*E;
end

if nargin > 3
    F = interp2(X, Y, F, Xq, Yq);
    if nargout > 1
        Fx = interp2(X, Y, Fx, Xq, Yq);
        Fy = interp2(X, Y, Fy, Xq, Yq);
    end
end

end
