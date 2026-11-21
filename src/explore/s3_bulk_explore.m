clear; root = piv_path('data');
res = piv_path('results');
load([res 'cyl_stats.mat'],'X','Y','D','mag','dt','U','V','stat_ok','uu','vv');
B = load([root 'Data/PIV_Data_Bulk.mat']); B = B.pivdata;
ub = B.dispX*mag/dt; vb = B.dispY*mag/dt;
bad = hypot(B.dispX,B.dispY)>8; ub(bad)=NaN; vb(bad)=NaN;
med = median(ub,3,'omitnan'); s = 1.4826*median(abs(ub-med),3,'omitnan'); o = abs(ub-med)>6*s; ub(o)=NaN; vb(o)=NaN;
Ub = mean(ub,3,'omitnan'); Vb = mean(vb,3,'omitnan'); urb = std(ub,0,3,'omitnan'); vrb = std(vb,0,3,'omitnan');
cols = 3:86;
prof = mean(Ub(:,cols),2,'omitnan'); yv = Y(:,1);
fprintf('Bulk: mean U over FOV %.4f m/s, std of profile %.4f; U at |y/D|<0.5: %.4f; top %.4f bottom %.4f\n', mean(prof), std(prof), mean(prof(abs(yv)<0.5)), prof(1), prof(end));
fprintf('Bulk: mean urms %.4f, vrms %.4f, Tu=%.2f%%\n', mean(mean(urb(:,cols))), mean(mean(vrb(:,cols))), 100*mean(mean(urb(:,cols)))/mean(prof));
fprintf('Bulk V mean %.4f m/s\n', mean(mean(Vb(:,cols),'omitnan')));
fprintf('profile (y/D : U):\n'); disp([yv(1:4:end) prof(1:4:end)]);
% convergence-ish: bulk time series at center
c = squeeze(ub(32,43,:)); fprintf('bulk point series std %.4f, lag-1 autocorr %.3f\n', std(c,'omitnan'), corr(c(1:end-1),c(2:end),'rows','complete'));
save([res 'bulk_stats.mat'],'Ub','Vb','urb','vrb','prof','yv');
% time-series exploration in wake
load([res 'cyl_velocity.mat'],'v','u');
[~,ix] = min(abs(X(1,:)-2.5)); [~,iy] = min(abs(Y(:,1)-0));
s1 = squeeze(v(iy,ix,:)); fprintf('wake probe x/D=%.2f y/D=%.2f valid %.2f\n', X(1,ix), Y(iy,1), mean(~isnan(s1)));
s1f = fillmissing(s1,'linear'); s1f = s1f-mean(s1f);
for L=1:6, fprintf('lag %d autocorr %.3f\n', L, corr(s1f(1:end-L),s1f(1+L:end))); end
[p,f] = pwelch(s1f, hann(512), 256, 512, 1);
[pm,j] = max(p(2:end)); fprintf('spectral peak at f = %.4f cycles/sample (period %.2f samples), peak/median = %.1f\n', f(j+1), 1/f(j+1), pm/median(p));
[~,ord] = sort(p,'descend'); fprintf('top peaks: %s\n', mat2str(round(f(ord(1:6))',4)));
fprintf('first 30 samples v: %s\n', mat2str(round(s1f(1:30)',3)));
