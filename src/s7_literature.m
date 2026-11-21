% Step 7: comparison with Tremblay (2001, Re=3900 DNS/LES + experiments) and
% analyses that the reference does not contain (triple decomposition, vortex tracking,
% shedding intermittency, convergence of Lr, top/bottom asymmetry, harmonic check).
clear; res = piv_path('results');
S = load([res 'cyl_stats.mat']); A = load([res 'analysis.mat']); load([res 'cyl_velocity.mat'],'u','v');
L = lit_tremblay();
U0 = A.U0; D = S.D; nu = S.nu;
fl = @(F) flipud(F);                         % rows -> y increasing
X = fl(S.X); Y = fl(S.Y); ok = fl(S.stat_ok);
U = fl(S.U); V = fl(S.V); uu = fl(S.uu); vv = fl(S.vv); uv = fl(S.uv);
u = flipud(u); v = flipud(v); nt = size(u,3);
xv = X(1,:); yv = Y(:,1);
hx = mean(diff(xv))*D; hy = mean(diff(yv))*D;
R.U0 = U0; R.xv = xv; R.yv = yv;

%% ---------- 1. centreline & profiles ----------
cl = @(F) interp1(yv,F,0);                    % value on y = 0
R.Ucl = cl(U)/U0; R.urms_cl = sqrt(cl(uu))/U0; R.vrms_cl = sqrt(cl(vv))/U0;
dn = xv>0.5 & ~isnan(R.Ucl);
xs = xv(dn); us = R.Ucl(dn);
iz = find(us(1:end-1)<0 & us(2:end)>=0,1,'last');
R.Lr_c = xs(iz) - us(iz)*(xs(iz+1)-xs(iz))/(us(iz+1)-us(iz)); R.Lr = R.Lr_c - 0.5;   % from rear surface
[R.Umin, j] = min(us); R.xUmin = xs(j);
[R.vrms_max, j] = max(R.vrms_cl(dn)); R.x_vmax = xs(j);
[R.urms_max_cl, j] = max(R.urms_cl(dn)); R.x_umax_cl = xs(j);
fprintf('Lr from rear = %.3f D (centre %.3f), Umin = %.3f at x/D = %.2f, vrms max %.3f at %.2f\n', R.Lr, R.Lr_c, R.Umin, R.xUmin, R.vrms_max, R.x_vmax);
% literature Lr from digitised centreline curves
lz = @(c) interp1(c(c(:,1)<2.5 & c(:,1)>1.2,2), c(c(:,1)<2.5 & c(:,1)>1.2,1), 0) - 0.5;
R.Lr_DNS = L.tab.Lr(3); R.Lr_LS = lz(L.cl.LS);
fprintf('digitised Lr (rear): LS %.2f; DNS table %.2f\n', R.Lr_LS, R.Lr_DNS);
% stations: equal x/D and equal x'/Lr (x' = x - D/2)
R.st_x = L.pu.x(1:3);
R.st_xLr = 0.5 + L.xLr*R.Lr;
prof = @(F,x) interp1(xv,F',x)';             % column at x
for q = 1:3
  for k = 1:2
    if k==1, x = R.st_x(q); else, x = R.st_xLr(q); end
    P.U{q,k} = prof(U,x)/U0; P.uu{q,k} = prof(uu,x)/U0^2; P.vv{q,k} = prof(vv,x)/U0^2; P.uv{q,k} = prof(uv,x)/U0^2;
  end
end
R.P = P;
% shear-layer peaks (upper & lower) at the equal-x'/Lr stations
for q=1:3
  p = P.uu{q,2}; up_ = yv>0.2 & yv<1.0; lo = yv<-0.2 & yv>-1.0;
  [R.uu_pk(q,1), j] = max(p.*up_); R.uu_pky(q,1) = yv(j);
  pl = p; pl(~lo) = -Inf; [R.uu_pk(q,2), j] = max(pl); R.uu_pky(q,2) = yv(j);
  p = P.uv{q,2}; pu_ = p; pu_(~up_) = Inf; [R.uv_pk(q,1), j] = min(pu_); R.uv_pky(q,1)=yv(j);
  pl = p; pl(~lo) = -Inf; [R.uv_pk(q,2), j] = max(pl); R.uv_pky(q,2)=yv(j);
end
disp('uu peaks [upper lower]:'); disp(R.uu_pk); disp('uv peaks [upper(min) lower(max)]:'); disp(R.uv_pk);

%% ---------- 2. triple decomposition (coherent / stochastic) ----------
bin = discretize(A.phi, linspace(-pi,pi,A.nbin+1)); use = A.use;
Uph = flipud(A.Uph); Vph = flipud(A.Vph); nph = A.nph;
su = zeros(size(U)); sv = su; suv = su; nn = su;           % stochastic sums
tu = su; tv = su; tuv = su; tn = su; mu_ = su; mv_ = su;   % totals over the used subset
for k = 1:nt
  if ~use(k), continue; end
  b = bin(k); uk = u(:,:,k); vk = v(:,:,k); g = ~isnan(uk) & ~isnan(vk) & ok;
  du = uk - Uph(:,:,b); dv = vk - Vph(:,:,b); g = g & ~isnan(du);
  du(~g)=0; dv(~g)=0; uk(~g)=0; vk(~g)=0;
  su = su + du.^2; sv = sv + dv.^2; suv = suv + du.*dv; nn = nn + g;
  mu_ = mu_ + uk; mv_ = mv_ + vk; tu = tu + uk.^2; tv = tv + vk.^2; tuv = tuv + uk.*vk;
end
mu_ = mu_./nn; mv_ = mv_./nn;
T.uu = tu./nn - mu_.^2; T.vv = tv./nn - mv_.^2; T.uv = tuv./nn - mu_.*mv_;        % total (subset)
Sx.uu = su./nn; Sx.vv = sv./nn; Sx.uv = suv./nn;                                 % stochastic
Cx.uu = T.uu - Sx.uu; Cx.vv = T.vv - Sx.vv; Cx.uv = T.uv - Sx.uv;                 % coherent
bad = nn < 200 | ~ok;
for F = {'uu','vv','uv'}, T.(F{1})(bad)=NaN; Sx.(F{1})(bad)=NaN; Cx.(F{1})(bad)=NaN; end
R.T = T; R.Sx = Sx; R.Cx = Cx; R.frac_c = (Cx.uu+Cx.vv)./(T.uu+T.vv);
R.cl_T = structfun(@(F) cl(F)/U0^2, T, 'UniformOutput', false);
R.cl_S = structfun(@(F) cl(F)/U0^2, Sx, 'UniformOutput', false);
R.cl_C = structfun(@(F) cl(F)/U0^2, Cx, 'UniformOutput', false);
fprintf('used snapshots %d of %d; coherent fraction of 2C TKE on axis at x/D=1.5,2,2.5: %s\n', nnz(use), nt, ...
  mat2str(round(100*interp1(xv, cl(R.frac_c), [1.5 2 2.5])),2));
% shear layer
[~,iy] = min(abs(yv-0.5)); R.frac_sl = R.frac_c(iy,:);

%% ---------- 3. vortex identification & tracking in phase-averaged fields ----------
th_ = linspace(0,2*pi,100); dA = hx*hy; trk = [];
for b = 1:A.nbin
  Ub = Uph(:,:,b); Vb = Vph(:,:,b); Ub(~ok)=NaN; Vb(~ok)=NaN;
  [ux,uy] = gradient(Ub,hx,hy); [vx,vy] = gradient(Vb,hx,hy);
  w = vx - uy; dsc = ((ux+vy)/2).^2 - (ux.*vy - uy.*vx);
  lci = sqrt(max(-dsc,0)); lci(isnan(lci)) = 0;
  thr = 0.15*max(lci(X>0.7)); reg = lci > thr & X > 0.6 & ok;
  CC = bwconncomp(reg,4);
  for r = 1:CC.NumObjects
    id = CC.PixelIdxList{r}; if numel(id) < 6, continue; end
    ww = w(id); sg = sign(sum(ww,'omitnan'));
    wt = abs(ww); wt(isnan(wt)) = 0;
    xc = sum(X(id).*wt)/sum(wt); yc = sum(Y(id).*wt)/sum(wt);
    G = sum(ww,'omitnan')*dA;                      % circulation of the swirling region [m^2/s]
    trk = [trk; b, sg, xc, yc, G/(U0*D), numel(id)]; %#ok<AGROW>
  end
end
R.trk = trk;   % [phase bin, sign, x/D, y/D, Gamma/(U0 D), npts]
% follow the strongest core of each sign through the phases (phase unwrapping)
phc = -pi + ((1:A.nbin)-0.5)*2*pi/A.nbin;
nb = A.nbin;
for sg = [1 -1]
  tt = trk(trk(:,2)==sg & trk(:,3)>0.8,:);
  xb = nan(nb,1); yb_ = xb; gb = xb;
  for b = 1:nb                                   % strongest core of this sign per phase
    q = tt(tt(:,1)==b,:); if isempty(q), continue; end
    [~,j] = max(abs(q(:,5))); xb(b) = q(j,3); yb_(b) = q(j,4); gb(b) = q(j,5);
  end
  best = [];
  for dirn = [1 -1]                              % phase may run either way round
    [~,b0] = min(xb); ord = mod(b0-1 + dirn*(0:nb-1), nb) + 1;
    xs_ = xb(ord); keep = 1;
    for k = 2:nb, if ~isnan(xs_(k)) && xs_(k) > xs_(keep(end)) - 0.05, keep(end+1) = k; else, break; end, end %#ok<AGROW>
    if numel(keep) > size(best,1)
      phu = (0:nb-1)'*2*pi/nb; phu = phu(keep); best = [phu xs_(keep) yb_(ord(keep)) gb(ord(keep))];
    end
  end
  c = polyfit(best(:,1), best(:,2), 1);
  key = sprintf('s%d', (sg>0)+1);
  R.track.(key) = best;                          % [phase since formation, x, y, Gamma]
  R.lambda.(key) = 2*pi*c(1);                    % streamwise advance per cycle (D)
  R.trackfit.(key) = c;
end
R.lambda_mean = mean([R.lambda.s1 R.lambda.s2]);
% cross-check: modes 1+2 form a travelling wave; spatial phase of v-component on the axis
m1 = flipud(A.modes(:,:,2,1)); m2 = flipud(A.modes(:,:,2,2));
th_cl = atan2(cl(m2), cl(m1)); gg = ~isnan(th_cl) & xv>0.9 & xv<2.6;
thu = unwrap(th_cl(gg)); c = polyfit(xv(gg), thu, 1); R.lambda_pod = 2*pi/abs(c(1));
R.th_cl = nan(size(xv)); R.th_cl(gg) = thu;
fprintf('wavelength from POD-mode phase on axis: %.2f D\n', R.lambda_pod);
fprintf('vortex wavelength lambda/D: + %.2f, - %.2f\n', R.lambda.s2, R.lambda.s1);
% with the shedding frequency implied by the (guessed) 1 Hz camera rate
fsh_guess = 0.4781;   % Hz: f/f_s = 0.478 (phase step, sign fixed by vortex travel, s9) with f_s = 1 Hz, the user's unverified estimate
R.Uc_U0 = R.lambda_mean*D*fsh_guess/U0; R.Uc_U0_pod = R.lambda_pod*D*fsh_guess/U0;
fprintf('POD-based Uc/U0 = %.2f\n', R.Uc_U0_pod);
fprintf('convection speed Uc/U0 = %.2f (uses f_shed = 0.478 Hz from the guessed 1 Hz rate)\n', R.Uc_U0);

%% ---------- 4. shedding intermittency (amplitude of POD modes 1-2) ----------
amp = A.amp; R.amp = amp;
[R.amp_pdf, R.amp_edges] = histcounts(amp, 0:0.05:3, 'Normalization','pdf');
mlag = 60; ac = zeros(mlag+1,1); a0 = amp - mean(amp);
for k = 0:mlag, ac(k+1) = sum(a0(1:end-k).*a0(1+k:end))/sum(a0.^2); end
R.amp_ac = ac; R.cyc_per_sample = 0.4781;   % f/f_s from the phase step (s9); alias order n = 0 assumed
q3 = quantile(amp,[1/3 2/3]);
grp = {amp<q3(1), amp>=q3(1) & amp<q3(2), amp>=q3(2)};
% centreline time series (interpolated to y=0)
[~,i0] = min(abs(yv)); if yv(i0)>0, i1=i0-1; else, i1=i0+1; end
wgt = abs(yv(i1))/(abs(yv(i0))+abs(yv(i1)));
ucl_t = squeeze(wgt*u(i0,:,:) + (1-wgt)*u(i1,:,:));        % nx x nt
for gq = 1:3
  m = mean(ucl_t(:,grp{gq}),2,'omitnan')/U0; R.cond.Ucl(:,gq) = m;
  s = m(dn); j = find(s(1:end-1)<0 & s(2:end)>=0,1,'last');
  if isempty(j), R.cond.Lr(gq) = NaN; else, R.cond.Lr(gq) = xs(j) - s(j)*(xs(j+1)-xs(j))/(s(j+1)-s(j)) - 0.5; end
  % profile at x/D = 1.06
  [~,ix] = min(abs(xv-1.06));
  R.cond.prof(:,gq) = mean(squeeze(u(:,ix,grp{gq})),2,'omitnan')/U0;
  R.cond.amp(gq) = mean(amp(grp{gq}));
end
fprintf('conditional Lr (low/mid/high amplitude): %s\n', mat2str(round(R.cond.Lr,3)));
% low-amplitude episodes: runs of consecutive samples with amp < 0.5
low = amp < 0.5; d = diff([0; low; 0]); st_ = find(d==1); en_ = find(d==-1)-1;
R.run_len = en_-st_+1; R.frac_low = mean(low);
fprintf('fraction amp<0.5: %.1f%%, mean run %.2f samples, max run %d\n', 100*R.frac_low, mean(R.run_len), max(R.run_len));
% expected run statistics if samples were independent with the same p
R.run_indep = 1/(1-R.frac_low);

%% ---------- 5. convergence of Lr with averaging length ----------
Nb = [10 20 40 80 160 320 640 1280];
for q = 1:numel(Nb)
  nbk = floor(nt/Nb(q)); lr = nan(nbk,1);
  for b = 1:nbk
    m = mean(ucl_t(:,(b-1)*Nb(q)+(1:Nb(q))),2,'omitnan')/U0; s = m(dn);
    j = find(s(1:end-1)<0 & s(2:end)>=0,1,'last');
    if ~isempty(j), lr(b) = xs(j) - s(j)*(xs(j+1)-xs(j))/(s(j+1)-s(j)) - 0.5; end
  end
  R.conv.N(q) = Nb(q); R.conv.mean(q) = mean(lr,'omitnan'); R.conv.std(q) = std(lr,'omitnan');
  R.conv.p05(q) = prctile(lr,5); R.conv.p95(q) = prctile(lr,95); R.conv.nb(q) = nbk;
end
disp([R.conv.N; R.conv.mean; R.conv.std]);

%% ---------- 6. top / bottom asymmetry ----------
B = A.bprof; yb = B.y; [yb,o] = sort(yb); Ubp = B.U(o);
g = ~isnan(Ubp) & abs(yb) < 0.8; c = polyfit(yb(g), Ubp(g), 1);
R.shear_beta = c(1)/U0;                 % (dU/dy) D / U0 of the approach flow at cylinder height
R.bulk_y = yb; R.bulk_U = Ubp/U0;
fprintf('approach-flow shear parameter beta = dU/dy D/U0 = %.3f\n', R.shear_beta);
% wake centre: y of minimum U between |y|<0.8 at each x
for i = 1:numel(xv)
  p = U(:,i); s = abs(yv)<0.8 & ~isnan(p);
  if nnz(s) < 5 || xv(i) < 0.6, R.yw(i) = NaN; continue; end
  [~,j] = min(p(s)); ys_ = yv(s); R.yw(i) = ys_(j);
end
% mean-V on axis and vortex circulation, upper vs lower
R.Vcl = cl(V)/U0;
R.G_up = min(trk(trk(:,2)<0,5)); R.G_lo = max(trk(trk(:,2)>0,5));
fprintf('peak circulation of tracked cores: clockwise %.3f, counter-clockwise %.3f (Gamma/U0D)\n', R.G_up, R.G_lo);

%% ---------- 7. harmonic content: aliased 2 f_shed in u on the axis ----------
nfft = 1024; w = 0.5-0.5*cos(2*pi*(0:nfft-1)'/(nfft-1)); R.sp.f = (0:nfft/2)'/nfft;
px = [1.25 1.5 2.0 2.5];
for q = 1:numel(px)
  [~,ix] = min(abs(xv-px(q)));
  su_ = fillmissing(ucl_t(ix,:)','linear','EndValues','nearest'); su_ = su_-mean(su_);
  [~,iy] = min(abs(yv-0)); sv_ = fillmissing(squeeze(v(iy,ix,:)),'linear','EndValues','nearest'); sv_ = sv_-mean(sv_);
  R.sp.Pu(:,q) = welch_(su_,nfft,w)/U0^2; R.sp.Pv(:,q) = welch_(sv_,nfft,w)/U0^2;
end
R.sp.x = px;
fpk = @(P,lo,hi) R.sp.f(find(P==max(P(R.sp.f>lo & R.sp.f<hi)),1));
fprintf('u-spectrum peaks near 0: %s; v-spectrum peaks near 0.5: %s\n', ...
  mat2str(arrayfun(@(q) fpk(R.sp.Pu(:,q),0.005,0.1),1:4),3), mat2str(arrayfun(@(q) fpk(R.sp.Pv(:,q),0.4,0.5),1:4),3));

save([res 'literature.mat'],'-struct','R');
disp('saved literature.mat');

function P = welch_(s,nfft,w)
  n = numel(s); step = nfft/2; nseg = floor((n-nfft)/step)+1; P = zeros(nfft/2+1,1);
  for k = 1:nseg
    seg = s((k-1)*step+(1:nfft)); seg = seg(:).*w; F = fft(seg); P = P + abs(F(1:nfft/2+1)).^2;
  end
  P = P/(nseg*sum(w.^2)); P(2:end-1) = 2*P(2:end-1);
end
