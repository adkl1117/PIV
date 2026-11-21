% Step 2: displacement -> velocity, outlier filtering, masks, statistics
clear; root = piv_path('data');
res = piv_path('results');
% ---- calibration (ring-peak circle fit on calibration image) ----
D_mm = 32.15; Dpx = 376.9; xc_px = 386.9; yc_row = 588.5;
mag = D_mm/Dpx*1e-3;                 % m/px
D = D_mm*1e-3;
nu = 9.17e-7; rho = 997.3;
P = load([root 'Data/PIV_Data_Cyl.mat']); P = P.pivdata;
dt = P.dt; cam = P.cam;
yc_up = cam(1) - yc_row;             % pixY axis points up
X = (P.pixX - xc_px)*mag/D;  Y = (P.pixY - yc_up)*mag/D;   % x/D, y/D, origin at cylinder centre
dx = P.dispX; dy = P.dispY; [ny,nx,nt] = size(dx);
nan0 = isnan(dx) | isnan(dy);
% sanity: mask centroid vs cylinder centre
nf = mean(nan0,3); m = nf>0.95 & hypot(X,Y)<0.6;
fprintf('NaN-core centroid: x/D=%.3f y/D=%.3f (should be ~0)\n', mean(X(m)), mean(Y(m)));
% ---- masks ----
R = hypot(X,Y);
mask_cyl = R <= 0.5 + 0.5*16*mag/D;                       % cylinder + half vector spacing
% ---- outlier detection ----
% (a) global: |disp| > 8 px (quarter of 32 px window) or beyond 6 robust-sigma of local time series
mag_d = hypot(dx,dy);
bad = mag_d > 8;
med_t = median(dx,3,'omitnan'); mad_t = 1.4826*median(abs(dx-med_t),3,'omitnan');
med_tv = median(dy,3,'omitnan'); mad_tv = 1.4826*median(abs(dy-med_tv),3,'omitnan');
bad = bad | abs(dx-med_t) > 6*max(mad_t,0.1) | abs(dy-med_tv) > 6*max(mad_tv,0.1);
nb1 = nnz(bad & ~nan0);
dx(bad)=NaN; dy(bad)=NaN;
% (b) universal outlier detection (Westerweel & Scarano 2005), 3x3, thr 2, eps 0.1 px
thr = 2; ep = 0.1; nb2 = 0;
for k = 1:nt
  for comp = 1:2
    if comp==1, F = dx(:,:,k); else, F = dy(:,:,k); end
    Fp = nan(ny+2,nx+2); Fp(2:end-1,2:end-1) = F;
    Nb = nan(ny,nx,8); c=0;
    for a=-1:1, for b=-1:1, if a==0&&b==0, continue; end; c=c+1; Nb(:,:,c)=Fp((2:end-1)+a,(2:end-1)+b); end, end
    Fm = median(Nb,3,'omitnan'); rm = median(abs(Nb-Fm),3,'omitnan');
    nval = sum(~isnan(Nb),3);
    r0 = abs(F-Fm)./(rm+ep); b0 = (r0>thr & nval>=3) | nval<2;
    if comp==1, B = b0; else, B = B | b0; end
  end
  B = B & ~isnan(dx(:,:,k));
  nb2 = nb2 + nnz(B);
  t1 = dx(:,:,k); t1(B)=NaN; dx(:,:,k)=t1; t2 = dy(:,:,k); t2(B)=NaN; dy(:,:,k)=t2;
end
dx(repmat(mask_cyl,1,1,nt)) = NaN; dy(repmat(mask_cyl,1,1,nt)) = NaN;
valid = ~isnan(dx);
fprintf('Vectors: %d total, NaN from processing %.1f%%, removed global/temporal %.2f%%, UOD %.2f%%, final valid %.1f%%\n', ...
   numel(dx), 100*mean(nan0(:)), 100*nb1/numel(dx), 100*nb2/numel(dx), 100*mean(valid(:)));
% ---- velocities ----
u = dx*mag/dt; v = dy*mag/dt;       % m/s
vf = mean(valid,3);                  % valid fraction per point
stat_ok = vf >= 0.5 & ~mask_cyl;     % enough samples for statistics
% sign check: upstream of cylinder, above centre -> v>0
V0 = mean(v,3,'omitnan'); up = X>-1.2 & X<-0.7;
fprintf('Mean v upstream above/below centre: %.4f / %.4f m/s\n', mean(V0(up & Y>0.3 & Y<0.8),'omitnan'), mean(V0(up & Y<-0.3 & Y>-0.8),'omitnan'));
% ---- statistics ----
U = mean(u,3,'omitnan'); V = mean(v,3,'omitnan');
up_ = u-U; vp_ = v-V;
uu = mean(up_.^2,3,'omitnan'); vv = mean(vp_.^2,3,'omitnan'); uv = mean(up_.*vp_,3,'omitnan');
Su = mean(up_.^3,3,'omitnan')./uu.^1.5; Sv = mean(vp_.^3,3,'omitnan')./vv.^1.5;
for F = {'U','V','uu','vv','uv','Su','Sv'}, eval([F{1} '(~stat_ok)=NaN;']); end
% free stream: upstream column band, outside |y|>1 (FOV left part above cylinder)
fs = X>-1 & X<-0.6 & abs(Y)>1.0 & stat_ok;
Uinf = mean(U(fs),'omitnan');
fsTop = X>1 & Y>1.3 & stat_ok; Utop = mean(U(fsTop),'omitnan');
ReD = Uinf*D/nu;
fprintf('U_inf (upstream, |y|>D) = %.4f m/s, U top region = %.4f m/s, Re_D = %.0f\n', Uinf, Utop, ReD);
% grid spacing (uniform)
hx = abs(X(1,2)-X(1,1))*D; hy = abs(Y(2,1)-Y(1,1))*D; sy = sign(Y(2,1)-Y(1,1));
ddx = @(F) gradient(F, hx, 1);  ddy = @(F) gradient(F, 1, hy*sy);
[dUdx,~] = gradient(U,hx,hy*sy); [~,dUdy] = gradient(U,hx,hy*sy);
[dVdx,~] = gradient(V,hx,hy*sy); [~,dVdy] = gradient(V,hx,hy*sy);
Wz = dVdx - dUdy;                                  % 1/s
Prod = -(uu.*dUdx + uv.*(dUdy+dVdx) + vv.*dVdy);  % 2D TKE production
k2 = 0.5*(uu+vv); k3 = 0.75*(uu+vv);
save([res 'cyl_stats.mat'],'X','Y','U','V','uu','vv','uv','Su','Sv','Wz','Prod','k2','k3','vf','stat_ok','mask_cyl', ...
     'Uinf','Utop','ReD','D','mag','dt','nu','rho','hx','hy','sy','nb1','nb2','nan0','-v7.3');
save([res 'cyl_velocity.mat'],'X','Y','u','v','D','mag','dt','Uinf','-v7.3');
fprintf('saved. nt=%d grid=%dx%d hx=%.2f mm\n', nt, ny, nx, hx*1e3);
