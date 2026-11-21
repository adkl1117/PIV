% Step 4: bulk reference, derived quantities, POD, phase averaging, spectra
clear; res = piv_path('results');
root = piv_path('data');
S = load([res 'cyl_stats.mat']); load([res 'cyl_velocity.mat'],'u','v');
X = S.X; Y = S.Y; D = S.D; mag = S.mag; dt = S.dt; nu = S.nu; rho = S.rho;
% ---------------- bulk (empty flume) ----------------
B = load([root 'Data/PIV_Data_Bulk.mat']); B = B.pivdata;
ub = B.dispX*mag/dt; vb = B.dispY*mag/dt; bad = hypot(B.dispX,B.dispY)>8;
med = median(ub,3,'omitnan'); s = 1.4826*median(abs(ub-med),3,'omitnan'); bad = bad | abs(ub-med)>6*s;
ub(bad)=NaN; vb(bad)=NaN; vfb = mean(~isnan(ub),3);
Ub = mean(ub,3,'omitnan'); Vb = mean(vb,3,'omitnan'); urb = std(ub,0,3,'omitnan'); vrb = std(vb,0,3,'omitnan');
uvb = mean((ub-Ub).*(vb-Vb),3,'omitnan');
okb = vfb>0.5; Ub(~okb)=NaN; Vb(~okb)=NaN; urb(~okb)=NaN; vrb(~okb)=NaN; uvb(~okb)=NaN;
bprof.y = Y(:,1); bprof.U = mean(Ub,2,'omitnan'); bprof.V = mean(Vb,2,'omitnan');
bprof.urms = mean(urb,2,'omitnan'); bprof.vrms = mean(vrb,2,'omitnan'); bprof.uv = mean(uvb,2,'omitnan');
bprof.Ustd_x = std(Ub,0,2,'omitnan'); nb_t = size(ub,3);
U0 = mean(bprof.U(abs(bprof.y)<=0.5));        % approach velocity at cylinder height
Ufov = mean(bprof.U,'omitnan');
fprintf('Bulk: U0(|y|<D/2)=%.4f, U_FOV=%.4f, Tu=%.2f%%, Tv=%.2f%%\n', U0, Ufov, ...
    100*mean(bprof.urms,'omitnan')/Ufov, 100*mean(bprof.vrms,'omitnan')/Ufov);
% ---------------- nondimensional numbers ----------------
mu = rho*nu; dp = 50e-6; rp = 1016; g = 9.81;
tau_p = rp*dp^2/(18*mu); ws = (rp-rho)*g*dp^2/(18*mu);
eps_ = U0^3/D; tau_eta = sqrt(nu/eps_); eta_K = (nu^3/eps_)^0.25; v_K = (nu*eps_)^0.25;
St_k = tau_p/tau_eta; St_D = tau_p*U0/D; ReD = U0*D/nu; Rep = ws*dp/nu;
nd = struct('U0',U0,'Ufov',Ufov,'ReD',ReD,'ReD_Utop',S.Utop*D/nu,'mu',mu,'tau_p',tau_p,'ws',ws, ...
    'eps',eps_,'tau_eta',tau_eta,'eta_K',eta_K,'v_K',v_K,'St_k',St_k,'St_D',St_D,'Rep',Rep,'Utop',S.Utop);
disp(nd)
% ---------------- derived mean quantities ----------------
ok = S.stat_ok; yv = Y(:,1); xv = X(1,:);
Ucl = interp1(yv, S.U, 0); Vcl = interp1(yv,S.V,0);
urmscl = interp1(yv,sqrt(S.uu),0); vrmscl = interp1(yv,sqrt(S.vv),0);
dn = xv>0.5 & ~isnan(Ucl); xs = xv(dn); us = Ucl(dn);
iz = find(us(1:end-1)<0 & us(2:end)>=0, 1, 'last');
if ~isempty(iz), Lr = xs(iz) - us(iz)*(xs(iz+1)-xs(iz))/(us(iz+1)-us(iz)); else, Lr = NaN; end
[Umin, im] = min(us);
fprintf('Recirculation: Lr (from centre) = %.3f D, from rear surface = %.3f D; min Ucl = %.3f U0 at x/D=%.2f\n', Lr, Lr-0.5, Umin/U0, xs(im));
vr = vrmscl(dn); [vmax, iv] = max(vr); fprintf('vrms max on centreline %.3f U0 at x/D=%.2f\n', vmax/U0, xs(iv));
stations = [0.75 1 1.25 1.5 2 2.5];
wake = zeros(numel(stations),4);
for q = 1:numel(stations)
  [~,ix] = min(abs(xv-stations(q))); p = S.U(:,ix); g_ = ~isnan(p);
  Ue = max(p(yv>0 & g_)); ud = Ue - min(p(abs(yv)<0.6 & g_)); half = Ue - ud/2;
  yy = yv(g_ & yv>0); pp = p(g_ & yv>0); [yy,o]=sort(yy); pp=pp(o);
  j = find(pp(1:end-1)<half & pp(2:end)>=half,1,'first');
  if isempty(j), bw=NaN; else, bw = yy(j)+(half-pp(j))*(yy(j+1)-yy(j))/(pp(j+1)-pp(j)); end
  wake(q,:) = [xv(ix) Ue/U0 ud/U0 bw];
  fprintf('x/D=%.2f: Ue=%.3f U0, deficit=%.3f U0, upper half-width y_1/2=%.3f D\n', wake(q,:));
end
% separation estimate: reversal of tangential velocity on ring r = 0.66 D, upper side
th = linspace(0,pi,181); rr = 0.66; xr = -rr*cos(th); yr = rr*sin(th);
Ut = interp2(X(end:-1:1,:),Y(end:-1:1,:),S.U(end:-1:1,:),xr,yr).*sin(th) + interp2(X(end:-1:1,:),Y(end:-1:1,:),S.V(end:-1:1,:),xr,yr).*cos(th);
js = find(Ut(1:end-1)>0 & Ut(2:end)<=0,1,'first');
if ~isempty(js), sep = th(js)*180/pi; else, sep = NaN; end
fprintf('Upper-side reversal of tangential velocity at r=%.2fD: theta = %.1f deg from front stagnation\n', rr, sep);
% ---------------- POD (snapshot via SVD) ----------------
nt = size(u,3);
pod_pts = S.stat_ok; pod_pts(:,[1 2 end]) = false; np_ = nnz(pod_pts);
% gappy data: fill missing fluctuations by spatial (Laplace) interpolation inside the domain
Uref = S.U; Vref = S.V; A = zeros(2*np_, nt);
for k = 1:nt
  fu = u(:,:,k)-Uref; fv = v(:,:,k)-Vref; fu(~pod_pts)=0; fv(~pod_pts)=0;
  hu = isnan(fu) & pod_pts; hv = isnan(fv) & pod_pts;
  fu(isnan(fu))=0; fv(isnan(fv))=0;
  if any(hu(:)), fu = regionfill(fu,hu); end
  if any(hv(:)), fv = regionfill(fv,hv); end
  A(:,k) = [fu(pod_pts); fv(pod_pts)];
end
gapfrac = 1-mean(mean(~isnan(u(repmat(pod_pts,1,1,nt)))));
fprintf('POD domain gap fraction (interpolated): %.1f%%\n', 100*gapfrac);
[Phi,Sig,Vt] = svd(A,'econ'); clear A
lam = diag(Sig).^2/nt; en_all = lam/sum(lam); a = Vt*Sig; clear Vt
fprintf('POD: %d points, energy modes 1-6: %s %%, cum 10: %.3f\n', np_, mat2str(round(100*en_all(1:6)',2)), sum(en_all(1:10)));
nm = 10; modes = nan([size(X) 2 nm]);
for k=1:nm
  t = nan(size(X)); t(pod_pts) = Phi(1:np_,k); modes(:,:,1,k)=t;
  t = nan(size(X)); t(pod_pts) = Phi(np_+1:end,k); modes(:,:,2,k)=t;
end
a = a(:,1:50); lam50 = lam(1:50); clear Phi Sig
% ---------------- phase from modes 1,2 ----------------
phi = atan2(a(:,2)/sqrt(2*lam(2)), a(:,1)/sqrt(2*lam(1)));
amp = sqrt(a(:,1).^2/(2*lam(1)) + a(:,2).^2/(2*lam(2)));
nbin = 12; edges = linspace(-pi,pi,nbin+1); bin = discretize(phi,edges);
use = amp > 0.5;
Uph = nan([size(X) nbin]); Vph = Uph; Wph = Uph; nph = zeros(nbin,1);
hx = S.hx; hy = S.hy; sy = S.sy;
for b=1:nbin
  sel = bin==b & use; nph(b)=nnz(sel);
  Ub_ = mean(u(:,:,sel),3,'omitnan'); Vb_ = mean(v(:,:,sel),3,'omitnan'); Ub_(~ok)=NaN; Vb_(~ok)=NaN;
  Uph(:,:,b)=Ub_; Vph(:,:,b)=Vb_;
  [dvx,~] = gradient(Vb_,hx,hy*sy); [~,duy] = gradient(Ub_,hx,hy*sy); Wph(:,:,b) = dvx-duy;
end
fprintf('phase bin counts: %s (used %d of %d)\n', mat2str(nph'), nnz(use), nt);
% ---------------- spectra (cycles per sample: acquisition rate unknown) ----------------
probes = [1.0 0.0; 1.5 0.0; 2.0 0.0; 2.5 0.0; 1.5 0.5; 2.0 -0.5; -0.7 1.2];
nfft = 512; w = 0.5-0.5*cos(2*pi*(0:nfft-1)'/(nfft-1));
spec.f = (0:nfft/2)'/nfft; spec.xy = probes; spec.valid = zeros(size(probes,1),1);
spec.Pu = zeros(nfft/2+1,size(probes,1)); spec.Pv = spec.Pu;
for i=1:size(probes,1)
  [~,ix] = min(abs(xv-probes(i,1))); [~,iy] = min(abs(yv-probes(i,2)));
  for comp=1:2
    if comp==1, s_ = squeeze(u(iy,ix,:)); else, s_ = squeeze(v(iy,ix,:)); end
    spec.valid(i) = mean(~isnan(s_));
    s_ = fillmissing(s_,'linear','EndValues','nearest'); s_ = s_-mean(s_);
    P = welch_(s_,nfft,w);
    if comp==1, spec.Pu(:,i)=P; else, spec.Pv(:,i)=P; end
  end
  spec.xy(i,:) = [xv(ix) yv(iy)];
end
spec.Pa1 = welch_(a(:,1)-mean(a(:,1)),nfft,w); spec.Pa3 = welch_(a(:,3)-mean(a(:,3)),nfft,w);
[~,j] = max(spec.Pa1(3:end)); fa = spec.f(j+2);
fprintf('Apparent peak of a1: %.4f cycles/sample\n', fa);
fsh = 0.21*U0/D; fprintf('Expected shedding f (St=0.21): with U0 %.3f Hz; with Utop %.3f Hz\n', fsh, 0.21*S.Utop/D);
ac = zeros(8,1); for L=1:8, ac(L) = corr(a(1:end-L,1),a(1+L:end,1)); end
fprintf('a1 autocorr lags 1..8: %s\n', mat2str(round(ac',3)));
% phase increment between consecutive samples
dphi = angle(exp(1i*diff(phi))); conv.dphi = dphi(use(1:end-1)&use(2:end));
fprintf('Phase step per sample: median %.1f deg, circular-mean %.1f deg, std %.1f deg\n', median(conv.dphi)*180/pi, angle(mean(exp(1i*conv.dphi)))*180/pi, std(conv.dphi)*180/pi);
% ---------------- convergence + pdfs ----------------
[~,ix] = min(abs(xv-1.5)); [~,iy] = min(abs(yv-0.45)); s_ = squeeze(u(iy,ix,:)); g_ = ~isnan(s_); s0 = s_; s0(~g_)=0;
cs = cumsum(s0)./cumsum(g_); cr = sqrt(max(cumsum(s0.^2)./cumsum(g_) - cs.^2,0));
conv.mean = cs; conv.rms = cr; conv.pt = [xv(ix) yv(iy)];
s2 = squeeze(v(iy,ix,:)); conv.pdf_u = s_(g_); conv.pdf_v = s2(~isnan(s2));
[~,ix2] = min(abs(xv-1.5)); [~,iy2] = min(abs(yv-0));
conv.pdf_v_cl = squeeze(v(iy2,ix2,:)); conv.pdf_u_cl = squeeze(u(iy2,ix2,:)); conv.pt_cl = [xv(ix2) yv(iy2)];
% ---------------- instantaneous snapshots at 4 phases ----------------
[~,ordA] = sort(amp,'descend'); cand = ordA(1:400); snaps = zeros(4,1);
for q=1:4, target = -pi + (q-1)*pi/2; [~,j] = min(abs(angle(exp(1i*(phi(cand)-target))))); snaps(q) = cand(j); end
snap_u = u(:,:,snaps); snap_v = v(:,:,snaps);
save([res 'analysis.mat'],'nd','bprof','Ub','Vb','urb','vrb','uvb','okb','nb_t','U0','Lr','Ucl','Vcl','urmscl','vrmscl', ...
  'stations','wake','sep','modes','en_all','lam50','a','phi','amp','use','nbin','Uph','Vph','Wph','nph','spec','fa','fsh','ac', ...
  'conv','snaps','snap_u','snap_v','pod_pts','Ut','th','rr','-v7.3');
disp('done');

function P = welch_(s,nfft,w)
  n = numel(s); nseg = floor((n-nfft)/(nfft/2))+1; P = zeros(nfft/2+1,1);
  for k=1:nseg, seg = s((k-1)*nfft/2+(1:nfft)).*w; F = fft(seg); P = P + abs(F(1:nfft/2+1)).^2; end
  P = P/(nseg*sum(w.^2)); P(2:end-1) = 2*P(2:end-1);
end
