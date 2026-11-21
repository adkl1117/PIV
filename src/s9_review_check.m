% Step 9: independent checks of the reviewer's claims (PIV/Review/审稿意见_PIV圆柱尾迹.md)
clear; res = piv_path('results');
root = piv_path('data');
S = load([res 'cyl_stats.mat']); A = load([res 'analysis.mat']); R = load([res 'literature.mat']);
load([res 'cyl_velocity.mat'],'u','v');
U0 = A.U0; D = S.D; nu = S.nu;
X = flipud(S.X); Y = flipud(S.Y); ok = flipud(S.stat_ok); u = flipud(u); v = flipud(v); nt = size(u,3);
xv = X(1,:); yv = Y(:,1); hx = mean(diff(xv))*D;
out = struct();

%% 1. direction of travel -> sign of the phase step
trk = R.trk; nb = A.nbin;
for sg = [-1 1]
  t = trk(trk(:,2)==sg & trk(:,3)>0.8,:); xb = nan(nb,1);
  for b = 1:nb, q = t(t(:,1)==b,:); if ~isempty(q), [~,j] = max(abs(q(:,5))); xb(b) = q(j,3); end, end
  dx = diff(xb); fprintf('sign %+d: core x per bin %s\n', sg, mat2str(round(xb',2)));
  fprintf('   median dx per +1 bin = %.3f D  (>0: vortices move downstream as phi increases)\n', median(dx(abs(dx)<0.5),'omitnan'));
end
dphi = angle(exp(1i*diff(A.phi))); g = A.use(1:end-1) & A.use(2:end);
m = angle(mean(exp(1i*dphi(g))));
fprintf('circular-mean phase step per sample = %.2f deg (resultant length %.2f)\n', m*180/pi, abs(mean(exp(1i*dphi(g)))));
fprintf('  -> if phi increases in physical time: frac(f/fs) = %.4f\n', mod(m,2*pi)/(2*pi));
fprintf('  -> if phi decreases in physical time: frac(f/fs) = %.4f\n', mod(-m,2*pi)/(2*pi));

%% 2. Lr from short blocks: last crossing vs. first crossing downstream of the minimum vs smoothed
[~,i0] = min(abs(yv)); if yv(i0)>0, i1=i0-1; else, i1=i0+1; end
w0 = abs(yv(i1))/(abs(yv(i0))+abs(yv(i1)));
ucl = squeeze(w0*u(i0,:,:) + (1-w0)*u(i1,:,:));
dn = find(xv>0.5); xs = xv(dn);
Lr_last  = @(s) xcross(xs, s, 'last');
Lr_first = @(s) xcross(xs(xs>1.27), s(xs>1.27), 'first');
Nb = [10 20 40 80 160 320]; T = zeros(numel(Nb),7);
for q = 1:numel(Nb)
  nbk = floor(nt/Nb(q)); L1 = nan(nbk,1); L2 = L1; L3 = L1;
  for b = 1:nbk
    s = mean(ucl(dn,(b-1)*Nb(q)+(1:Nb(q))),2,'omitnan')/U0;
    L1(b) = Lr_last(s); s2 = s(xs>1.27); L2(b) = Lr_first(s);
    L3(b) = xcross(xs, smoothdata(s,'movmean',5,'omitnan'), 'last');
  end
  T(q,:) = [Nb(q) mean(L1,'omitnan') std(L1,'omitnan') mean(L2,'omitnan') std(L2,'omitnan') mean(L3,'omitnan') std(L3,'omitnan')];
end
Lfull = Lr_last(mean(ucl(dn,:),2,'omitnan')/U0);
disp('N | last: mean std | first-after-min: mean std | smoothed-last: mean std   (all from rear surface)');
disp(round([T(:,1) T(:,2:7)-[0.5 0 0.5 0 0.5 0]],3)); fprintf('full record Lr = %.3f\n', Lfull-0.5);

%% 3. valid fraction vs phase bin and vs local gradient
bin = discretize(A.phi, linspace(-pi,pi,nb+1));
pts = [1.5 0.45; 1.0 0.5; 2.0 0; 1.5 -0.45];
for p = 1:size(pts,1)
  [~,ix] = min(abs(xv-pts(p,1))); [~,iy] = min(abs(yv-pts(p,2)));
  vk = squeeze(~isnan(u(iy,ix,:))); vf = arrayfun(@(b) mean(vk(bin==b & A.use)), 1:nb);
  fprintf('valid fraction at (%.2f,%.2f): overall %.2f, per phase bin %.2f..%.2f\n', xv(ix), yv(iy), mean(vk), min(vf), max(vf));
end

%% 4. amplitude autocorrelation: is there a lag-0 discontinuity (noise floor)?
ac = R.amp_ac(1:8)'; fprintf('A autocorr lags 0..7: %s\n', mat2str(round(ac,3)));
c = polyfit(1:4, log(ac(2:5)), 1); fprintf('exp fit lags 1-4: rho(0) extrapolated = %.2f, ratio per lag = %.2f -> e-fold %.2f samples\n', exp(c(2)), exp(c(1)), -1/c(1));
a1 = A.a(:,1); r1 = arrayfun(@(L) corr(a1(1:end-L),a1(1+L:end)), 1:8); fprintf('a1 autocorr lags 1..8: %s\n', mat2str(round(r1,2)));

%% 5. dissipation: 2/3 law on the stochastic part, and integral scale
Uph = flipud(A.Uph); sel_y = find(yv>-0.55 & yv<0.75); sel_x = find(xv>1.0 & xv<2.6);
rr = 1:6; DLL = zeros(size(rr)); nn = DLL; Ruu = zeros(1,30); nR = Ruu; var_ = 0; nv = 0;
for k = 1:nt
  if ~A.use(k), continue; end
  f_ = u(sel_y,sel_x,k) - Uph(sel_y,sel_x,bin(k));
  for j = rr
    d = f_(:,1+j:end) - f_(:,1:end-j); g_ = ~isnan(d); DLL(j) = DLL(j) + sum(d(g_).^2); nn(j) = nn(j) + nnz(g_);
  end
  f0 = f_ - mean(f_(:),'omitnan');
  for j = 0:29
    p_ = f0(:,1:end-j).*f0(:,1+j:end); g_ = ~isnan(p_); Ruu(j+1) = Ruu(j+1) + sum(p_(g_)); nR(j+1) = nR(j+1) + nnz(g_);
  end
end
DLL = DLL./nn; Ruu = Ruu./nR; rho_ = Ruu/Ruu(1); z = find(rho_<=0,1); if isempty(z), z = numel(rho_); end
L11 = trapz((0:z-1)*hx, rho_(1:z)); urms = sqrt(Ruu(1));
eps23 = (DLL/2.0).^1.5./(rr*hx);
fprintf('u''''_rms = %.4f m/s (%.2f U0), L11 = %.1f mm (%.2f D)\n', urms, urms/U0, L11*1e3, L11/D);
fprintf('eps from 2/3 law at r = 1..6 h: %s m^2/s^3\n', mat2str(eps23,2));
fprintf('eps = 0.5 u^3/L11 = %.2e; eps = u^3/L11 = %.2e; worksheet u_b^3/D = %.2e\n', 0.5*urms^3/L11, urms^3/L11, 0.0542^3/D);   % u_b with b = 60 mm (s6)
for e = [2e-4 9e-4 5.0e-3], fprintf('  eps %.1e: eta = %.0f um, tau = %.1f ms\n', e, (nu^3/e)^0.25*1e6, sqrt(nu/e)*1e3); end

%% 6. measured bulk profile power law with bed at -60 mm, and Q with m in 3..10
y = A.bprof.y*D; ub = A.bprof.U; gq = ~isnan(ub); y = y(gq); ub = ub(gq); [y,o] = sort(y); ub = ub(o);
zz = y + 0.060; low = y < 0.02;                % lower part of the FOV, below the velocity maximum
cp = polyfit(log(zz(low)), log(ub(low)), 1); fprintf('profile power law z^(1/m), lower FOV: 1/m = %.2f (m = %.1f)\n', cp(1), 1/cp(1));
b = 0.060; h = 0.118;   % width 60 mm (user, 2026-11-21)
ybot = min(y); Lfov = max(y)-ybot; I_fov = trapz(y,ub); u_bot = ub(1); u_top = ub(end);
pt = polyfit(y(end-5:end), ub(end-5:end), 1);
rng(1); N = 2e5;
for mlo = [4 3 2]
  yb = -0.063 + 0.006*rand(N,1); gb = ybot - yb; gt = h - Lfov - gb; mb = mlo + (10-mlo)*rand(N,1);
  Q = b*(u_bot*gb.*mb./(mb+1) + I_fov + u_top*gt + rand(N,1).*0.5*pt(1).*gt.^2) .* (1+0.015*(2*rand(N,1)-1)) .* (1+0.01*(2*rand(N,1)-1)) .* ...
      (1 - 2*min((0.5+0.5*rand(N,1))*h, b/2)./(b*((5+4*rand(N,1))+1)));
  pr = prctile(Q,[2.5 50 97.5])*1e3; fprintf('m in [%d,10]: Q median %.2f l/s, 95%% %.2f..%.2f\n', mlo, pr(2), pr(1), pr(3));
end
fprintf('Q with m = 3 fixed (others central): ratio to m = 7: %.3f\n', (u_bot*0.0228*3/4 + I_fov + u_top*(h-Lfov-0.0228))/(u_bot*0.0228*7/8 + I_fov + u_top*(h-Lfov-0.0228)));

%% 7. is the 11 % inflow turbulence real?
B = load([root 'Data/PIV_Data_Bulk.mat']); B = B.pivdata; mag = S.mag; dt = S.dt;
ubk = B.dispX*mag/dt; ubk(hypot(B.dispX,B.dispY)>8) = NaN;
fm = squeeze(mean(mean(ubk(10:55,5:80,:),1,'omitnan'),2,'omitnan'));
fprintf('bulk: frame-averaged u rms = %.2f mm/s; local u rms = %.2f mm/s; 0.1 px noise = %.2f mm/s\n', std(fm)*1e3, mean(std(ubk(10:55,5:80,:),0,3,'omitnan'),'all','omitnan')*1e3, 0.1*mag/dt*1e3);
s1 = squeeze(ubk(32,20,:)); s2 = squeeze(ubk(32,60,:)); gg = ~isnan(s1)&~isnan(s2);
fprintf('bulk: correlation of u over %.0f mm = %.2f\n', 40*16*mag*1e3, corr(s1(gg),s2(gg)));
fprintf('Re with blockage-corrected U_c = 0.0745 (b = 60 mm): %.0f\n', 0.0745*D/nu);
save([res 'review_check.mat'],'T','Lfull','eps23','L11','urms');

function xc = xcross(x, s, which)
  x = x(:); s = s(:); j = find(s(1:end-1)<0 & s(2:end)>=0, 1, which);
  if isempty(j), xc = NaN; else, xc = x(j) - s(j)*(x(j+1)-x(j))/(s(j+1)-s(j)); end
end
