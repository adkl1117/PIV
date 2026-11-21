% Step 1: magnification from cylinder silhouette (D = 32.15 mm)
root = piv_path('data');
res  = piv_path('results');
D_mm = 32.15;
d = dir([root 'Samples/*.TIF']);
A = zeros(1040,1392);
for k=1:numel(d), A = A + double(imread([root 'Samples/' d(k).name])); end
A = A/numel(d);                                  % mean particle image (particles outside, none inside)
G = imgaussfilt(A,3);
c0 = [387 590];                                  % rough centre (x, row)
th = (0:1:359)'*pi/180; r = (100:0.5:280)';
rim = nan(size(th));
for i=1:numel(th)
  xs = c0(1)+r*cos(th(i)); ys = c0(2)-r*sin(th(i));
  ok = xs>1 & xs<1392 & ys>1 & ys<1040; if ~all(ok), continue; end
  p = interp2(G,xs,ys);
  pin = median(p(r<150)); pout = median(p(r>230));
  if pout-pin < 3, continue; end
  j = find(p > (pin+pout)/2, 1, 'first');       % first crossing of mid-level going outward
  if ~isempty(j), rim(i) = r(j); end
end
ok = ~isnan(rim);
for it=1:4   % robust algebraic circle fit with outlier rejection
  X = c0(1)+rim(ok).*cos(th(ok)); Y = c0(2)-rim(ok).*sin(th(ok));
  M = [X Y ones(size(X))]; b = -(X.^2+Y.^2); p = M\b;
  xc = -p(1)/2; yc = -p(2)/2; R = sqrt(xc^2+yc^2-p(3));
  e = sqrt((X-xc).^2+(Y-yc).^2)-R; s = 1.4826*median(abs(e));
  idx = find(ok); ok(idx(abs(e)>2.5*s)) = false;
end
X = c0(1)+rim(ok).*cos(th(ok)); Y = c0(2)-rim(ok).*sin(th(ok));
e = sqrt((X-xc).^2+(Y-yc).^2)-R;
% horizontal / vertical diameters separately (ellipse check)
Dh = (max(X)-min(X)); Dv = (max(Y)-min(Y));
mag = D_mm/(2*R);                                % mm per px
fprintf('Circle fit: xc=%.1f px, row=%.1f px, R=%.2f px, D=%.1f px, N=%d pts, rms=%.2f px\n', xc,yc,R,2*R,sum(ok),rms(e));
fprintf('Extent: Dh=%.1f px, Dv=%.1f px\n', Dh, Dv);
fprintf('Magnification = %.5f mm/px  (%.2f px/mm)\n', mag, 1/mag);
% sensitivity: +-2 px in D
fprintf('Uncertainty for +-3 px edge error: %.5f .. %.5f mm/px (+-%.1f%%)\n', D_mm/(2*R+3), D_mm/(2*R-3), 100*3/(2*R));
cal.mag_mm_per_px = mag; cal.xc_px = xc; cal.yc_row = yc; cal.R_px = R; cal.D_mm = D_mm; cal.rms = rms(e);
save([res 'calibration.mat'],'cal');
f = figure('Visible','off','Position',[0 0 1100 820]);
imagesc(G); colormap gray; axis image; hold on; caxis(prctile(G(:),[1 99.7]));
plot(X,Y,'c.','MarkerSize',6); t=linspace(0,2*pi,400); plot(xc+R*cos(t),yc+R*sin(t),'r-','LineWidth',1.2); plot(xc,yc,'r+','MarkerSize',14);
title(sprintf('Cylinder edge fit: D = %.1f px  \rightarrow  %.4f mm/px', 2*R, mag));
exportgraphics(f,[res 'fig00_calibration.png'],'Resolution',130);
