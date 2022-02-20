function gs = green_to_mat(js, gs)

gs.name = js.name; 
fprintf('%s:\n\tphi=%d\n\tnStg=%d\n', gs.name, gs.phi, gs.nStg);

if ~isfield(gs, 'Stg')
    t = linspace(0,1,gs.nStg)';
    gs.Stg = 0.02*250.^t;
end

Nx = 400;
x0 = max(0*js.Dj,js.X(1,1));
x1 = min(20*js.Dj,js.X(1,end));
gs.Xg = linspace(x0,x1,Nx)';

Nr = 150;
r0 = js.Y(1,1);
r1 = min(3.5*js.Dj,js.Y(end,1));
gs.Rg = linspace(r0,r1,Nr)';

[gs.X gs.R] = meshgrid(gs.Xg, gs.Rg);

gs = green(js, gs);

matfile = sprintf('GREEN/%s_green_%d',gs.name, gs.phi);
save(matfile, '-struct', 'gs');
fprintf('done\n');
end
