% exploration: bulk profile shape, feasible bed range, profile-law fits
clear; res = piv_path('results');
A = load([res 'analysis.mat'],'bprof','Ub'); S = load([res 'cyl_stats.mat'],'D');
D = S.D*1e3; h = 118; bp = A.bprof;
y = bp.y*D; u = bp.U; g = ~isnan(u); y = y(g); u = u(g); [y,o] = sort(y); u = u(o);
fprintf('FOV y range: %.1f .. %.1f mm (relative to cylinder centre)\n', min(y), max(y));
fprintf('feasible bed position y_b: %.1f .. %.1f mm  (FOV inside water column, depth %g mm)\n', max(y)-h, min(y), h);
[um,im] = max(u); fprintf('u_max = %.4f m/s at y = %.1f mm\n', um, y(im));
fprintf('profile (every 3rd):\n'); disp([y(1:3:end) u(1:3:end)]);
% power-law fit with free bed position on points below the maximum
sel = y < y(im)-3;
ys = y(sel); us = u(sel);
best = [];
for yb = linspace(max(y)-h, min(y)-0.5, 400)
  z = ys - yb; X = [ones(size(z)) log(z)]; c = X\log(us); r = log(us) - X*c;
  best = [best; yb, 1/c(2), exp(c(1)), sum(r.^2)]; %#ok<AGROW>
end
[~,j] = min(best(:,4)); fprintf('power law: best y_b = %.1f mm, 1/m exponent -> m = %.2f, rss=%.3g (rss at edges: %.3g / %.3g)\n', best(j,1), best(j,2), best(j,4), best(1,4), best(end,4));
% log law fit, free y_b
bestL = [];
for yb = linspace(max(y)-h, min(y)-0.5, 400)
  z = (ys - yb)*1e-3; X = [ones(size(z)) log(z)]; c = X\us; r = us - X*c;
  bestL = [bestL; yb, c(2)*0.41, exp(-c(1)/c(2)), sum(r.^2)]; %#ok<AGROW>
end
[~,j] = min(bestL(:,4)); fprintf('log law: best y_b = %.1f mm, u* = %.4f m/s, z0 = %.3g m, rss=%.3g (edges %.3g / %.3g)\n', bestL(j,1), bestL(j,2), bestL(j,3), bestL(j,4), bestL(1,4), bestL(end,4));
fprintf('dip hypothesis z_max/h = 0.80..0.90 -> y_b = %.1f .. %.1f mm\n', y(im)-0.9*h, y(im)-0.8*h);
fprintf('cylinder at mid-depth -> y_b = %.1f mm\n', -h/2);
% column-to-column scatter
U2 = A.Ub(:,3:85); s = std(U2,0,2,'omitnan'); fprintf('mean column-to-column std of profile: %.4f m/s\n', mean(s,'omitnan'));
