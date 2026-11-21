% Step 8: figures for the literature comparison (Tremblay 2001) and the new analyses
clear; close all;
res = piv_path('results'); fig = [res 'figures/'];
R = load([res 'literature.mat']); A = load([res 'analysis.mat'],'Wph','nbin','nph'); S = load([res 'cyl_stats.mat'],'X','Y','stat_ok');
L = lit_tremblay();
X = flipud(S.X); Y = flipud(S.Y); ok = flipud(S.stat_ok); xv = R.xv; yv = R.yv;
hex = @(h) sscanf(h(2:end),'%2x%2x%2x',[1 3])/255;
cat = {'#2a78d6','#eb6834','#1baf7a','#eda100','#e87ba4','#008300','#6250d6','#e34948'};
C = cellfun(hex,cat,'UniformOutput',false);
ink = hex('#0b0b0b'); ink2 = hex('#52514e'); mut = hex('#8a8984'); nodata = hex('#cfcdc6'); cylc = hex('#3a3936');
grid_c = [0.88 0.88 0.86];
seqk = {'#f4f8fd','#cde2fb','#9ec5f4','#6da7ec','#3987e5','#256abf','#184f95','#0d366b'};
seq = interp1(linspace(0,1,numel(seqk)), cell2mat(cellfun(hex,seqk','UniformOutput',false)), linspace(0,1,256));
divk = {'#0d366b','#256abf','#6da7ec','#cde2fb','#f0efec','#f8d0cf','#ee8a89','#d23c3b','#8a1f1f'};
div = interp1(linspace(0,1,numel(divk)), cell2mat(cellfun(hex,divk','UniformOutput',false)), linspace(0,1,256));
set(groot,'defaultAxesFontSize',10,'defaultTextColor',ink,'defaultAxesXColor',ink2,'defaultAxesYColor',ink2, ...
    'defaultAxesTickDir','out','defaultAxesBox','off','defaultFigureColor','w','defaultAxesLineWidth',0.8, ...
    'defaultLineLineWidth',2,'defaultAxesTitleFontWeight','normal','defaultAxesTitleFontSizeMultiplier',1.05, ...
    'defaultLegendBox','off','defaultLegendFontSize',9);
gr = @() set(gca,'XGrid','on','YGrid','on','GridColor',grid_c,'GridAlpha',1);
th = linspace(0,2*pi,200); cyl = @() fill(0.5*cos(th),0.5*sin(th),cylc,'EdgeColor','none');
mir = @(c) [-flipud(c(:,1)) flipud(c(:,2)); c];          % symmetric extension
mira = @(c) [-flipud(c(:,1)) -flipud(c(:,2)); c];        % antisymmetric extension (u'v')
ours = C{1}; dns = ink2; ls = C{2}; ow = C{3};
Lr = R.Lr; LrD = L.tab.Lr(3); LrLS = R.Lr_LS;

%% L01 centreline
f = figure('Position',[0 0 1500 430]); t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; hold on
plot(L.cl.DNS(:,1),L.cl.DNS(:,2),'-','Color',dns,'LineWidth',1.4);
plot(L.cl.LS(:,1),L.cl.LS(:,2),'+','Color',ls,'MarkerSize',8,'LineWidth',1.2);
plot(L.cl.OW(:,1),L.cl.OW(:,2),'x','Color',ow,'MarkerSize',9,'LineWidth',1.4);
g = xv>0.5; plot(xv(g),R.Ucl(g),'-','Color',ours,'LineWidth',2.4);
yline(0,'Color',mut,'LineWidth',0.8);
plot(0.5+Lr,0,'o','MarkerFaceColor',ours,'MarkerEdgeColor','w','MarkerSize',8);
plot(0.5+LrD,0,'o','MarkerFaceColor',dns,'MarkerEdgeColor','w','MarkerSize',7);
xlim([0.5 4.5]); ylim([-0.45 1]); gr(); xlabel('x/D'); ylabel('\langle u\rangle / U_0 on y = 0');
legend({'Tremblay DNS, Re 3900','Lourenco & Shih PIV','Ong & Wallace HW','this experiment, Re \approx 1900'},'Location','southeast');
title('(a) Centreline mean velocity');
text(0.5+Lr+0.05,-0.08,sprintf('L_r = %.2f D',Lr),'Color',ink,'FontSize',9);
nexttile; hold on
plot((L.cl.DNS(:,1)-0.5)/LrD,L.cl.DNS(:,2),'-','Color',dns,'LineWidth',1.4);
plot((L.cl.LS(:,1)-0.5)/LrLS,L.cl.LS(:,2),'+','Color',ls,'MarkerSize',8,'LineWidth',1.2);
plot((xv(g)-0.5)/Lr,R.Ucl(g),'-','Color',ours,'LineWidth',2.4); yline(0,'Color',mut,'LineWidth',0.8); xline(1,'Color',mut,'LineWidth',0.8);
xlim([0 2.6]); ylim([-0.45 1]); gr(); xlabel('x''/L_r   (x'' = x - D/2)'); ylabel('\langle u\rangle / U_0');
title('(b) Same, rescaled with each case''s L_r'); legend({'DNS','L&S','this exp.'},'Location','southeast');
nexttile; hold on
plot(xv(g),R.vrms_cl(g),'-','Color',ours,'LineWidth',2.4); plot(xv(g),R.urms_cl(g),'-','Color',C{4},'LineWidth',2.4);
xs4 = L.pu.x; vD = cellfun(@(c) sqrt(c(1,2)), L.vv.DNS); vL = cellfun(@(c) sqrt(c(1,2)), L.vv.LS);
uD = cellfun(@(c) sqrt(c(1,2)), L.uu.DNS); uL = cellfun(@(c) sqrt(c(1,2)), L.uu.LS);
plot(xs4,vD,'s','MarkerFaceColor',dns,'MarkerEdgeColor','w','MarkerSize',8);
plot(xs4,vL,'+','Color',ls,'MarkerSize',9,'LineWidth',1.4);
plot(xs4,uD,'d','MarkerFaceColor','w','MarkerEdgeColor',dns,'MarkerSize',7,'LineWidth',1.2);
plot(xs4,uL,'x','Color',ls,'MarkerSize',9,'LineWidth',1.4);
xline(0.5+Lr,':','Color',ours,'LineWidth',1.2); xline(0.5+LrD,':','Color',dns,'LineWidth',1.2);
xlim([0.5 3.2]); ylim([0 0.9]); gr(); xlabel('x/D'); ylabel('rms / U_0 on y = 0');
legend({'v_{rms} this exp.','u_{rms} this exp.','v_{rms} DNS','v_{rms} L&S','u_{rms} DNS','u_{rms} L&S'},'Location','northeast','NumColumns',2);
title('(c) Centreline fluctuations (dotted: bubble closure)');
exportgraphics(f,[fig 'L01_centreline_vs_lit.png'],'Resolution',140); close(f);

%% L02 mean-velocity profiles: equal x/D (top) and equal x'/Lr (bottom)
f = figure('Position',[0 0 1400 780]); t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for k = 1:2
  for q = 1:3
    nexttile; hold on
    cD = mir(L.pu.DNS{q}); cL = mir(L.pu.LS{q});
    plot(cD(:,1),cD(:,2),'-','Color',dns,'LineWidth',1.4); plot(cL(:,1),cL(:,2),'+','Color',ls,'MarkerSize',7,'LineWidth',1.1);
    plot(yv,R.P.U{q,k},'-','Color',ours,'LineWidth',2.4);
    xlim([-1.6 1.6]); ylim([-0.45 1.75]); gr(); yline(0,'Color',mut,'LineWidth',0.8);
    if k==1, ttl = sprintf('x/D = %.2f (all)', L.pu.x(q)); else
      ttl = sprintf('x''/L_r = %.2f: lit. x/D = %.2f, ours x/D = %.2f', L.xLr(q), L.pu.x(q), R.st_xLr(q)); end
    title(ttl,'FontSize',9.5); xlabel('y/D'); if q==1, ylabel('\langle u\rangle / U_0'); end
    if k==1 && q==1, legend({'Tremblay DNS','Lourenco & Shih','this exp.'},'Location','north'); end
  end
end
title(t,'Mean streamwise velocity profiles. Top row: same x/D.  Bottom row: same x''/L_r (Tremblay Fig. 4.26 scaling)','FontSize',11,'Color',ink);
exportgraphics(f,[fig 'L02_profiles_vs_lit.png'],'Resolution',140); close(f);

%% L03 Reynolds-stress profiles at equal x'/Lr
f = figure('Position',[0 0 1400 1000]); t = tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
nm = {'uu','vv','uv'}; lab = {'\langle u''u''\rangle / U_0^2','\langle v''v''\rangle / U_0^2','\langle u''v''\rangle / U_0^2'};
yl_ = {[0 0.32],[0 0.65],[-0.16 0.16]};
for r = 1:3
  for q = 1:3
    nexttile; hold on
    if r==3, cD = mira(L.(nm{r}).DNS{q}); cL = mira(L.(nm{r}).LS{q}); else, cD = mir(L.(nm{r}).DNS{q}); cL = mir(L.(nm{r}).LS{q}); end
    fill([cD(:,1); flipud(cD(:,1))],[cD(:,2)-0.02; flipud(cD(:,2)+0.02)],dns,'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
    plot(cD(:,1),cD(:,2),'-','Color',dns,'LineWidth',1.4); plot(cL(:,1),cL(:,2),'+','Color',ls,'MarkerSize',7,'LineWidth',1.1);
    plot(yv,R.P.(nm{r}){q,2},'-','Color',ours,'LineWidth',2.4);
    xlim([-1.6 1.6]); ylim(yl_{r}); gr(); if r==3, yline(0,'Color',mut,'LineWidth',0.8); end
    if r==1, title(sprintf('x''/L_r = %.2f  (ours x/D = %.2f, lit. %.2f)', L.xLr(q), R.st_xLr(q), L.pu.x(q)),'FontSize',9.5); end
    if q==1, ylabel(lab{r}); end; if r==3, xlabel('y/D'); end
    if r==1 && q==1, legend({'Tremblay DNS','Lourenco & Shih','this exp.'},'Location','northwest'); end
  end
end
title(t,'Reynolds stresses at equal x''/L_r. Grey band: ±0.02 digitisation uncertainty of the DNS curve','FontSize',11,'Color',ink);
exportgraphics(f,[fig 'L03_stresses_vs_lit.png'],'Resolution',140); close(f);

%% N01 triple decomposition
f = figure('Position',[0 0 1500 460]); t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
ax = nexttile; set(ax,'Color',nodata); hold on
F = 100*R.frac_c; F(~ok) = NaN; pcolor(X,Y,F); shading interp; colormap(ax,seq); caxis([0 100]); cb = colorbar; cb.Label.String = 'coherent share of k_{2C} [%]';
contour(X,Y,F,[50 50],'Color',C{2},'LineWidth',1.4); cyl(); axis image; xlim([min(xv) max(xv)]); ylim([min(yv) max(yv)]);
xlabel('x/D'); ylabel('y/D'); title('(a) Share of in-plane TKE carried by the shedding cycle');
nexttile; hold on; g = xv>0.6;
plot(xv(g),R.cl_T.vv(g),'-','Color',C{1},'LineWidth',2.4); plot(xv(g),R.cl_C.vv(g),'--','Color',C{1},'LineWidth',1.8); plot(xv(g),R.cl_S.vv(g),':','Color',C{1},'LineWidth',2);
plot(xv(g),R.cl_T.uu(g),'-','Color',C{4},'LineWidth',2.4); plot(xv(g),R.cl_C.uu(g),'--','Color',C{4},'LineWidth',1.8); plot(xv(g),R.cl_S.uu(g),':','Color',C{4},'LineWidth',2);
gr(); xlabel('x/D'); ylabel('stress / U_0^2 on y = 0'); xlim([0.6 max(xv)]);
legend({'v''v'' total','v''v'' coherent','v''v'' stochastic','u''u'' total','u''u'' coherent','u''u'' stochastic'},'Location','northwest','NumColumns',2);
title('(b) Axis: coherent + stochastic = total');
nexttile; hold on
rT = R.cl_T.vv./R.cl_T.uu; rS = R.cl_S.vv./R.cl_S.uu;
plot(xv(g),rT(g),'-','Color',C{1},'LineWidth',2.4); plot(xv(g),rS(g),'-','Color',C{3},'LineWidth',2.4);
rD = cellfun(@(c) c(1,2), L.vv.DNS)./cellfun(@(c) c(1,2), L.uu.DNS); rL = cellfun(@(c) c(1,2), L.vv.LS)./cellfun(@(c) c(1,2), L.uu.LS);
plot(L.pu.x,rD,'s','MarkerFaceColor',dns,'MarkerEdgeColor','w','MarkerSize',8); plot(L.pu.x,rL,'+','Color',ls,'MarkerSize',9,'LineWidth',1.4);
yline(1,'Color',mut,'LineWidth',0.8); gr(); xlim([0.6 3.2]); ylim([0 9]); xlabel('x/D'); ylabel('\langle v''v''\rangle / \langle u''u''\rangle on y = 0');
legend({'this exp. total','this exp. stochastic part','Tremblay DNS (total)','L&S (total)'},'Location','northwest');
title('(c) Axis anisotropy: total vs. stochastic');
exportgraphics(f,[fig 'N01_triple_decomposition.png'],'Resolution',140); close(f);

%% N02 vortex tracking
trk = R.trk; f = figure('Position',[0 0 1500 880]); t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile; set(ax,'Color',nodata); hold on; b0 = 4;
W = flipud(A.Wph(:,:,b0))*R.xv(1)*0; W = flipud(A.Wph(:,:,b0))*0.03215/R.U0; W(~ok)=NaN;
pcolor(X,Y,W); shading interp; colormap(ax,div); caxis([-4 4]); cb = colorbar; cb.Label.String = '\langle\omega_z\rangle_\phi D/U_0';
q = trk(trk(:,1)==b0,:);
for r = 1:size(q,1), plot(q(r,3),q(r,4),'o','MarkerSize',6+18*abs(q(r,5))/max(abs(trk(:,5))),'MarkerEdgeColor',ink,'LineWidth',1.4); end
cyl(); axis image; xlim([min(xv) max(xv)]); ylim([min(yv) max(yv)]); xlabel('x/D'); ylabel('y/D');
title(sprintf('(a) Phase bin %d: detected cores (circle size ~ |\\Gamma|)', b0));
nexttile; hold on
pp = trk(trk(:,2)<0,:); scatter(pp(:,3),pp(:,4),20+200*abs(pp(:,5))/max(abs(trk(:,5))),C{1},'filled','MarkerEdgeColor','w');
pp = trk(trk(:,2)>0,:); scatter(pp(:,3),pp(:,4),20+200*abs(pp(:,5))/max(abs(trk(:,5))),C{8},'filled','MarkerEdgeColor','w');
cyl(); axis image; xlim([min(xv) max(xv)]); ylim([min(yv) max(yv)]); gr(); xlabel('x/D'); ylabel('y/D');
legend({'clockwise (\omega<0, from upper side)','counter-clockwise (\omega>0, from lower side)'},'Location','southoutside','NumColumns',2);
title('(b) Core positions over all 12 phases (vortex-street paths)');
nexttile; hold on
s1 = R.track.s1; s2 = R.track.s2;
plot(s1(:,1)/(2*pi),s1(:,2),'o','MarkerFaceColor',C{1},'MarkerEdgeColor','w','MarkerSize',8);
plot(s2(:,1)/(2*pi),s2(:,2),'o','MarkerFaceColor',C{8},'MarkerEdgeColor','w','MarkerSize',8);
for s = {s1,s2}, z = s{1}; c = polyfit(z(:,1)/(2*pi),z(:,2),1); xx = linspace(0,max(z(:,1))/(2*pi),10); plot(xx,polyval(c,xx),'-','Color',ink2,'LineWidth',1.2); end
gr(); xlabel('phase since the core appears / 2\pi  [shedding cycles]'); ylabel('core x/D');
legend({'clockwise core','counter-clockwise core','linear fit'},'Location','northwest');
title(sprintf('(c) Core advance per cycle: %.2f D (cw), %.2f D (ccw); POD-mode phase: %.2f D', R.lambda.s1, R.lambda.s2, R.lambda_pod),'FontSize',9.5);
nexttile; hold on
pp = trk(trk(:,2)<0,:); plot(pp(:,3),abs(pp(:,5)),'o','MarkerFaceColor',C{1},'MarkerEdgeColor','w','MarkerSize',7);
pp = trk(trk(:,2)>0,:); plot(pp(:,3),abs(pp(:,5)),'o','MarkerFaceColor',C{8},'MarkerEdgeColor','w','MarkerSize',7);
gr(); xlabel('core x/D'); ylabel('|\Gamma| / (U_0 D) of swirling region'); xline(0.5+Lr,':','Color',ink2,'LineWidth',1.2);
legend({'clockwise','counter-clockwise','bubble closure'},'Location','northwest');
title('(d) Circulation of phase-averaged cores vs. position');
exportgraphics(f,[fig 'N02_vortex_tracking.png'],'Resolution',140); close(f);

%% N03 intermittency of shedding
f = figure('Position',[0 0 1500 860]); t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on; n0 = 1:400; cps = R.cyc_per_sample; cps = 0.478;
plot(n0*cps, R.amp(n0),'-','Color',C{1},'LineWidth',1.2); yline(0.5,'--','Color',C{2},'LineWidth',1.2);
gr(); xlabel('time [shedding cycles]  (sample number \times 0.48, assumes n = 0)'); ylabel('shedding amplitude A = |a_1+ia_2| (norm.)');
title('(a) Amplitude of the shedding mode, first 400 samples'); legend({'A(t)','A = 0.5 (excluded from phase average)'},'Location','northeast');
nexttile; hold on; e = R.amp_edges; c_ = (e(1:end-1)+e(2:end))/2;
bar(c_,R.amp_pdf,1,'FaceColor',C{1},'EdgeColor','w','FaceAlpha',0.9);
aa = linspace(0,3,200); plot(aa,2*aa.*exp(-aa.^2),'-','Color',C{2},'LineWidth',2);
gr(); xlabel('A'); ylabel('pdf'); legend({'measured','Rayleigh (random narrow-band noise)'},'Location','northeast');
title(sprintf('(b) Amplitude PDF; A < 0.5 in %.1f %% of samples', 100*R.frac_low));
nexttile; hold on; lag = (0:numel(R.amp_ac)-1)*cps;
plot(lag,R.amp_ac,'-o','Color',C{1},'MarkerFaceColor',C{1},'MarkerSize',3,'LineWidth',1.6); yline(0,'Color',mut,'LineWidth',0.8);
yline(exp(-1),':','Color',ink2,'LineWidth',1); gr(); xlabel('lag [shedding cycles]  (assumes n = 0)'); ylabel('autocorrelation of A');
title('(c) Memory of the shedding amplitude (envelope modulation)');
nexttile; hold on; g = xv>0.5 & xv<2.55;
plot(L.cl.DNS(:,1),L.cl.DNS(:,2),'-','Color',dns,'LineWidth',1.4); plot(L.cl.LS(:,1),L.cl.LS(:,2),'+','Color',ls,'MarkerSize',7,'LineWidth',1.1);
plot(xv(g),R.cond.Ucl(g,1),'-','Color',C{3},'LineWidth',2.2); plot(xv(g),R.cond.Ucl(g,2),'-','Color',C{1},'LineWidth',2.2); plot(xv(g),R.cond.Ucl(g,3),'-','Color',C{7},'LineWidth',2.2);
yline(0,'Color',mut,'LineWidth',0.8); xlim([0.5 3.5]); ylim([-0.45 1]); gr(); xlabel('x/D'); ylabel('\langle u\rangle / U_0 on y = 0');
legend({'Tremblay DNS','Lourenco & Shih',sprintf('weak shedding (A<%.2f): L_r = %.2f D',quantile(R.amp,1/3),R.cond.Lr(1)), ...
  sprintf('middle third: L_r = %.2f D',R.cond.Lr(2)), sprintf('strong shedding: L_r = %.2f D',R.cond.Lr(3))},'Location','southeast');
title('(d) Centreline velocity conditioned on shedding strength');
exportgraphics(f,[fig 'N03_shedding_intermittency.png'],'Resolution',140); close(f);

%% N04 convergence of Lr (revised after review: two estimators, snapshot sampling, n = 0 assumed)
V = load([res 'review_check.mat']); T = V.T; ncy = T(:,1)'*R.cyc_per_sample;
f = figure('Position',[0 0 960 500]); hold on
fill([ncy fliplr(ncy)],[T(:,2)'-0.5-T(:,3)' fliplr(T(:,2)'-0.5+T(:,3)')],C{1},'FaceAlpha',0.15,'EdgeColor','none');
fill([ncy fliplr(ncy)],[T(:,4)'-0.5-T(:,5)' fliplr(T(:,4)'-0.5+T(:,5)')],C{2},'FaceAlpha',0.15,'EdgeColor','none');
plot(ncy,T(:,2)-0.5,'-o','Color',C{1},'MarkerFaceColor',C{1},'LineWidth',2);
plot(ncy,T(:,4)-0.5,'-s','Color',C{2},'MarkerFaceColor',C{2},'LineWidth',2);
yline(Lr,'--','Color',ink2,'LineWidth',1.2);
set(gca,'XScale','log'); gr(); xlabel('block length [shedding cycles]  (sample number \times 0.48, assumes n = 0)'); ylabel('L_r / D  (from rear surface)');
legend({'\pm1 std, last zero crossing','\pm1 std, first crossing after the minimum','block mean, last zero crossing','block mean, first crossing after the minimum',sprintf('full record: %.2f D',Lr)},'Location','northeast');
title('L_r from blocks of consecutive snapshots: the short-block bias depends on the estimator');
exportgraphics(f,[fig 'N04_Lr_convergence.png'],'Resolution',140); close(f);

%% N05 asymmetry from the free surface / bed
f = figure('Position',[0 0 1500 450]); t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; hold on; yb = R.bulk_y; plot(R.bulk_U,yb,'-','Color',C{1},'LineWidth',2.4);
fill([0.6 1.3 1.3 0.6],[-0.5 -0.5 0.5 0.5],cylc,'FaceAlpha',0.15,'EdgeColor','none');
gr(); xlabel('empty-flume \langle u\rangle / U_0'); ylabel('y/D'); ylim([min(yb) max(yb)]);
title(sprintf('(a) Approach flow is sheared: \\beta = (dU/dy)D/U_0 = %.3f', R.shear_beta),'FontSize',9.5);
text(0.62,0,' cylinder height','Color',ink,'FontSize',9);
nexttile; hold on; g = xv>0.6;
tmp = load([res 'cyl_stats.mat'],'uv'); UV = flipud(tmp.uv); yw = nan(size(xv));
for i = find(g)
  p = UV(:,i); s = abs(yv) < 0.5 & ~isnan(p); ys_ = yv(s); ps_ = p(s);
  j = find(ps_(1:end-1) > 0 & ps_(2:end) <= 0, 1);          % +uv below, -uv above
  if ~isempty(j), yw(i) = ys_(j) - ps_(j)*(ys_(j+1)-ys_(j))/(ps_(j+1)-ps_(j)); end
end
yw = smoothdata(yw,'movmean',5,'omitnan'); yw(xv<1.3) = NaN;   % inside the bubble uv changes sign several times
plot(xv(g),yw(g),'-','Color',C{1},'LineWidth',2.4); plot(xv(g),R.Vcl(g),'-','Color',C{2},'LineWidth',2.4);
yline(0,'Color',mut,'LineWidth',0.8); gr(); xlabel('x/D'); ylabel('y_w/D   or   \langle v\rangle/U_0');
legend({'wake centre y_w (zero of \langle u''v''\rangle)','\langle v\rangle / U_0 on y = 0'},'Location','best');
save([res 'wake_centre.mat'],'yw');
title('(b) Wake centre and mean cross-flow on the axis');
nexttile; hold on
vals = [R.uu_pk(:,1)'; R.uu_pk(:,2)']; vals2 = [abs(R.uv_pk(:,1))'; abs(R.uv_pk(:,2))'];
bb = bar(1:3,[vals' vals2'],'grouped','EdgeColor','w'); cols = {C{1},hex('#9ec5f4'),C{2},hex('#f3b48f')};
for k=1:4, bb(k).FaceColor = cols{k}; end
set(gca,'XTick',1:3,'XTickLabel',arrayfun(@(z) sprintf('x''/L_r=%.2f',z),L.xLr,'UniformOutput',false)); gr();
ylabel('peak value / U_0^2'); legend({'u''u'' upper','u''u'' lower','|u''v''| upper','|u''v''| lower'},'Location','northeast','NumColumns',2);
title('(c) Upper vs. lower shear-layer peaks');
exportgraphics(f,[fig 'N05_asymmetry.png'],'Resolution',140); close(f);

%% N06 harmonic content of the phase-averaged velocity on the axis (independent of f_s)
P = load([res 'analysis.mat'],'Uph','Vph'); Uph = flipud(P.Uph); Vph = flipud(P.Vph);
nb = size(Uph,3); phc = -pi + ((1:nb)-0.5)*2*pi/nb;
ua = squeeze(interp1(yv,Uph,0))/R.U0; va = squeeze(interp1(yv,Vph,0))/R.U0;   % nx x nb
H = @(F,k) 2*abs(F*exp(-1i*k*phc(:)))/nb;                                    % harmonic amplitude
Hu1 = H(ua,1); Hu2 = H(ua,2); Hv1 = H(va,1); Hv2 = H(va,2);
f = figure('Position',[0 0 1800 450]); t = tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
nexttile; hold on; [~,ix] = min(abs(xv-2.0));
pp = [phc phc(1)+2*pi]*180/pi; uu_ = ua(ix,:)-mean(ua(ix,:)); vv_ = va(ix,:)-mean(va(ix,:));
plot(pp,[vv_ vv_(1)],'-o','Color',C{1},'MarkerFaceColor',C{1},'MarkerSize',5,'LineWidth',2);
plot(pp,[uu_ uu_(1)],'-o','Color',C{4},'MarkerFaceColor',C{4},'MarkerSize',5,'LineWidth',2);
yline(0,'Color',mut,'LineWidth',0.8); gr(); xlim([-180 180]); set(gca,'XTick',-180:90:180);
xlabel('shedding phase \phi [deg]'); ylabel('phase-averaged fluctuation / U_0'); legend({'\langle v\rangle_\phi - \langle v\rangle','\langle u\rangle_\phi - \langle u\rangle'},'Location','southwest');
title(sprintf('(a) Axis, x/D = %.2f: v: one clean cycle; u: mixed f and 2f', xv(ix)));
nexttile; hold on; Hu1 = Hu1(:)'; Hu2 = Hu2(:)'; Hv1 = Hv1(:)'; Hv2 = Hv2(:)'; g = xv>0.7 & ~isnan(Hu1);
plot(xv(g),Hv1(g),'-','Color',C{1},'LineWidth',2.4); plot(xv(g),Hv2(g),'--','Color',C{1},'LineWidth',1.8);
plot(xv(g),Hu1(g),'-','Color',C{4},'LineWidth',2.4); plot(xv(g),Hu2(g),'--','Color',C{4},'LineWidth',1.8);
gr(); xlabel('x/D on y = 0'); ylabel('harmonic amplitude / U_0'); xlim([0.7 max(xv)]);
legend({'v: 1st harmonic (f_{shed})','v: 2nd harmonic','u: 1st harmonic','u: 2nd harmonic (2 f_{shed})'},'Location','northwest');
title('(b) Harmonic amplitudes along the axis');
nexttile; hold on
plot(xv(g),Hu2(g)./Hu1(g),'-','Color',C{4},'LineWidth',2.4); plot(xv(g),Hv2(g)./Hv1(g),'-','Color',C{1},'LineWidth',2.4);
yline(1,'Color',mut,'LineWidth',0.8); gr(); set(gca,'YScale','log'); xlabel('x/D on y = 0'); ylabel('2nd / 1st harmonic'); xlim([0.7 max(xv)]);
legend({'u','v'},'Location','northeast'); title('(c) On y = 0, u keeps a strong f_{shed} part');
nexttile; hold on
E = exp(-1i*reshape(phc,1,1,[])); Hu1m = 2*abs(sum(Uph.*E,3))/nb/R.U0;
ymin = nan(size(xv)); gg = find(xv>1.2);
for i = gg, p = Hu1m(:,i); s = abs(yv)<0.6 & ~isnan(p); ys_ = yv(s); if nnz(s) < 5, continue; end; [~,j] = min(p(s)); ymin(i) = ys_(j); end
ymin = smoothdata(ymin,'movmean',5,'omitnan'); W2 = load([res 'wake_centre.mat']);
plot(xv(gg),ymin(gg),'-','Color',C{4},'LineWidth',2.4); plot(xv(gg),W2.yw(gg),'-','Color',C{1},'LineWidth',2.4);
yline(0,'Color',mut,'LineWidth',0.8); gr(); xlabel('x/D'); ylabel('y/D'); ylim([-0.3 0.3]); xlim([1.2 max(xv)]);
legend({'y where u has no f_{shed} part (dynamic centre)','zero of \langle u''v''\rangle (mean centre)'},'Location','southwest');
title('(d) Where a symmetric wake would put the axis');
R2.ymin = ymin;
exportgraphics(f,[fig 'N06_harmonics_axis.png'],'Resolution',140); close(f);
R2.Hu1 = Hu1; R2.Hu2 = Hu2; R2.Hv1 = Hv1; R2.Hv2 = Hv2; save([res 'harmonics.mat'],'-struct','R2');
fprintf('axis x=1.5/2/2.5: u 2nd/1st = %s, v 2nd/1st = %s\n', mat2str(round(interp1(xv,Hu2./Hu1,[1.5 2 2.5]),2)), mat2str(round(interp1(xv,Hv2./Hv1,[1.5 2 2.5]),2)));
fprintf('dynamic centre y at x=1.5/2/2.5: %s; mean centre: %s\n', mat2str(round(interp1(xv,R2.ymin,[1.5 2 2.5]),3)), mat2str(round(interp1(xv,W2.yw,[1.5 2 2.5]),3)));
fprintf('u 1st/2nd harmonic amplitude on axis at x=2: %.3f / %.3f; v 1st %.3f\n', interp1(xv,Hu1,2), interp1(xv,Hu2,2), interp1(xv,Hv1,2));
disp('figures done');
