% Step 5: all figures for the cylinder-wake PIV analysis
clear; close all;
res = piv_path('results');
fig = [res 'figures/']; if ~exist(fig,'dir'), mkdir(fig); end
root = piv_path('data');
S = load([res 'cyl_stats.mat']); A = load([res 'analysis.mat']);
U0 = A.U0; D = S.D;
% ---- flip to Y increasing (meshgrid form) ----
fl = @(F) flipud(F);
X = fl(S.X); Y = fl(S.Y);
U = fl(S.U)/U0; V = fl(S.V)/U0; uu = fl(S.uu)/U0^2; vv = fl(S.vv)/U0^2; uv = fl(S.uv)/U0^2;
W = fl(S.Wz)*D/U0; k3 = fl(S.k3)/U0^2; Pr = fl(S.Prod)*D/U0^3; Su = fl(S.Su); Sv = fl(S.Sv); vf = fl(S.vf);
ok = fl(S.stat_ok);
xl = [min(X(:)) max(X(:))]; yl = [min(Y(:)) max(Y(:))];
% ---- palette (reference palette, light mode) ----
hex = @(h) sscanf(h(2:end),'%2x%2x%2x',[1 3])/255;
cat = {'#2a78d6','#eb6834','#1baf7a','#eda100','#e87ba4','#008300','#6250d6','#e34948'};
C = cellfun(hex,cat,'UniformOutput',false);
ink = hex('#0b0b0b'); ink2 = hex('#52514e'); nodata = hex('#cfcdc6'); cylc = hex('#3a3936');
seqk = {'#f4f8fd','#cde2fb','#9ec5f4','#6da7ec','#3987e5','#256abf','#184f95','#0d366b'};
seq = interp1(linspace(0,1,numel(seqk)), cell2mat(cellfun(hex,seqk','UniformOutput',false)), linspace(0,1,256));
divk = {'#0d366b','#256abf','#6da7ec','#cde2fb','#f0efec','#f8d0cf','#ee8a89','#d23c3b','#8a1f1f'};
div = interp1(linspace(0,1,numel(divk)), cell2mat(cellfun(hex,divk','UniformOutput',false)), linspace(0,1,256));
set(groot,'defaultAxesFontSize',10,'defaultTextColor',ink,'defaultAxesXColor',ink2,'defaultAxesYColor',ink2, ...
    'defaultAxesTickDir','out','defaultAxesBox','off','defaultFigureColor','w','defaultAxesLineWidth',0.8, ...
    'defaultLineLineWidth',2,'defaultAxesTitleFontWeight','normal','defaultAxesTitleFontSizeMultiplier',1.1);
th = linspace(0,2*pi,200);
cyl = @() fill(0.5*cos(th),0.5*sin(th),cylc,'EdgeColor','none');
    function_holder = []; %#ok<NASGU>

%% F01: raw image + valid fraction
I = double(imread([root 'Samples/cyl_000000.T000.D000.P000.H000.LA.TIF']));
f = figure('Position',[0 0 1500 560]);
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; mag = S.mag; xc = 386.9; yr = 588.5;
xi = ((1:1392)-xc)*mag/D; yi = ((1040-(1:1040))-(1040-yr))*mag/D;
imagesc(xi,yi,I); set(gca,'YDir','normal'); colormap(gca,gray); caxis(prctile(I(:),[1 99.6])); axis image; hold on
plot(0.5*cos(th),0.5*sin(th),'Color',C{2},'LineWidth',1.2);
plot(X(1:4:end,1:4:end),Y(1:4:end,1:4:end),'.','Color',C{1},'MarkerSize',4);
xlabel('x/D'); ylabel('y/D'); title('Raw particle image (frame A, #0) with every 4th vector position');
nexttile; pcolor(X,Y,100*vf); shading flat; axis image; hold on; cyl(); colormap(gca,seq); cb=colorbar; cb.Label.String='valid vectors [%]';
caxis([0 100]); xlabel('x/D'); ylabel('y/D'); title(sprintf('Valid-vector fraction after filtering (mean %.1f %%)', 100*mean(vf(:))));
contour(X,Y,vf,[0.5 0.5],'Color',C{2},'LineWidth',1.2);
text(-0.1,-0.85,'laser shadow','Color',ink,'FontSize',9,'HorizontalAlignment','center');
exportgraphics(f,[fig 'F01_raw_valid.png'],'Resolution',140); close(f);

%% helper for field plots
function_field = @(ax,F,cl,cm,ttl,cbl) []; %#ok<NASGU>

%% F02: mean velocity magnitude + streamlines
Um = hypot(U,V); Um(~ok) = NaN;
f = figure('Position',[0 0 1100 760]); ax = axes; set(ax,'Color',nodata); hold on
pcolor(X,Y,Um); shading interp; colormap(seq); caxis([0 1.7]); cb = colorbar; cb.Label.String = '|\langle U \rangle| / U_0';
U2 = U; V2 = V; U2(isnan(U2))=0; V2(isnan(V2))=0;
hs = streamslice(X,Y,U2,V2,2.2); set(hs,'Color',[ink 0.75],'LineWidth',0.7);
contour(X,Y,U,[0 0],'--','Color',C{2},'LineWidth',1.6);
cyl(); axis image; xlim(xl); ylim(yl); xlabel('x/D'); ylabel('y/D');
title(sprintf('Time-averaged velocity magnitude and streamlines (N = %d, Re_D \\approx %.0f)', size(A.a,1), A.nd.ReD));
text(A.Lr,0.08,sprintf(' L_r = %.2f D',A.Lr),'Color',C{2},'FontSize',10,'FontWeight','bold');
plot(A.Lr,0,'o','MarkerFaceColor',C{2},'MarkerEdgeColor','w','MarkerSize',8);
text(xl(1)+0.05,yl(1)+0.07,'gray = no valid data (cylinder, laser shadow)','Color',ink2,'FontSize',8.5);
exportgraphics(f,[fig 'F02_mean_streamlines.png'],'Resolution',140); close(f);

%% F03: mean U, V
f = figure('Position',[0 0 1500 560]); tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; set(gca,'Color',nodata); hold on; pcolor(X,Y,U); shading interp; colormap(gca,div); caxis([-1.7 1.7]); cb=colorbar; cb.Label.String='\langle u \rangle / U_0';
contour(X,Y,U,[0 0],'k-','LineWidth',1.2); contour(X,Y,U,[0.5 1 1.5],':','Color',ink2,'LineWidth',0.8);
cyl(); axis image; xlim(xl); ylim(yl); xlabel('x/D'); ylabel('y/D'); title('Streamwise mean velocity (black: \langle u \rangle = 0)');
nexttile; set(gca,'Color',nodata); hold on; pcolor(X,Y,V); shading interp; colormap(gca,div); caxis([-0.5 0.5]); cb=colorbar; cb.Label.String='\langle v \rangle / U_0';
contour(X,Y,U,[0 0],'k-','LineWidth',1.2); cyl(); axis image; xlim(xl); ylim(yl); xlabel('x/D'); title('Transverse mean velocity');
exportgraphics(f,[fig 'F03_mean_U_V.png'],'Resolution',140); close(f);

%% F04: mean vorticity + vectors
f = figure('Position',[0 0 1100 760]); ax=axes; set(ax,'Color',nodata); hold on
Wm = W; Wm(~ok)=NaN; pcolor(X,Y,Wm); shading interp; colormap(div); caxis([-6 6]); cb=colorbar; cb.Label.String='\langle\omega_z\rangle D / U_0';
s = 3; quiver(X(1:s:end,1:s:end),Y(1:s:end,1:s:end),U(1:s:end,1:s:end),V(1:s:end,1:s:end),1.3,'Color',[ink 0.8],'LineWidth',0.6);
cyl(); axis image; xlim(xl); ylim(yl); xlabel('x/D'); ylabel('y/D');
title('Mean out-of-plane vorticity with mean velocity vectors (red: CCW, blue: CW)');
exportgraphics(f,[fig 'F04_mean_vorticity.png'],'Resolution',140); close(f);

%% F05: Reynolds stresses
f = figure('Position',[0 0 1500 1050]); tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
flds = {uu, vv, uv, k3}; cls = {[0 0.4],[0 0.6],[-0.15 0.15],[0 0.6]}; cms = {seq,seq,div,seq};
lbl = {'\langle u''u''\rangle / U_0^2','\langle v''v''\rangle / U_0^2','\langle u''v''\rangle / U_0^2','k / U_0^2,  k = 3/4(\langle u''^2\rangle+\langle v''^2\rangle)'};
ttl = {'Streamwise normal stress','Transverse normal stress','Reynolds shear stress','Turbulent kinetic energy (3D estimate)'};
for q=1:4
  nexttile; set(gca,'Color',nodata); hold on; F = flds{q}; F(~ok)=NaN; pcolor(X,Y,F); shading interp; colormap(gca,cms{q}); caxis(cls{q});
  cb=colorbar; cb.Label.String=lbl{q}; contour(X,Y,U,[0 0],'-','Color',ink2,'LineWidth',0.9); cyl(); axis image; xlim(xl); ylim(yl);
  xlabel('x/D'); ylabel('y/D'); title(ttl{q});
end
exportgraphics(f,[fig 'F05_reynolds_stresses.png'],'Resolution',140); close(f);

%% F06: production + skewness
f = figure('Position',[0 0 1500 520]); tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; set(gca,'Color',nodata); hold on; F=Pr; F(~ok)=NaN; pcolor(X,Y,F); shading interp; colormap(gca,div); caxis([-1.2 1.2]);
cb=colorbar; cb.Label.String='P_k D / U_0^3'; cyl(); axis image; xlim(xl); ylim(yl); xlabel('x/D'); ylabel('y/D'); title('TKE production (in-plane terms)');
nexttile; set(gca,'Color',nodata); hold on; F=Su; F(~ok)=NaN; pcolor(X,Y,F); shading interp; colormap(gca,div); caxis([-1.5 1.5]);
cb=colorbar; cb.Label.String='S_u'; cyl(); axis image; xlim(xl); ylim(yl); xlabel('x/D'); title('Skewness of u''');
nexttile; set(gca,'Color',nodata); hold on; F=Sv; F(~ok)=NaN; pcolor(X,Y,F); shading interp; colormap(gca,div); caxis([-1.5 1.5]);
cb=colorbar; cb.Label.String='S_v'; cyl(); axis image; xlim(xl); ylim(yl); xlabel('x/D'); title('Skewness of v''');
exportgraphics(f,[fig 'F06_production_skewness.png'],'Resolution',140); close(f);

%% F07: centreline
xv = S.X(1,:); f = figure('Position',[0 0 1100 520]); hold on
dn = xv>0.5;
yline(0,'-','Color',[0.8 0.8 0.78],'LineWidth',1,'HandleVisibility','off');
p1=plot(xv(dn),A.Ucl(dn)/U0,'-','Color',C{1}); p2=plot(xv(dn),A.urmscl(dn)/U0,'-','Color',C{2}); p3=plot(xv(dn),A.vrmscl(dn)/U0,'-','Color',C{3});
p4=plot(xv(dn),A.Vcl(dn)/U0,'-','Color',C{4});
xline(A.Lr,'--','Color',ink2,'LineWidth',1,'HandleVisibility','off');
text(A.Lr+0.03,1.15,sprintf('L_r = %.2f D',A.Lr),'Color',ink2);
xlabel('x/D  (y/D = 0)'); ylabel('velocity / U_0'); xlim([0.5 xl(2)]);
legend([p1 p2 p3 p4],{'\langle u\rangle','u_{rms}','v_{rms}','\langle v\rangle'},'Location','northwest','Box','off');
title('Wake centreline: mean velocity and fluctuation intensity'); grid on; set(gca,'GridColor',[0.85 0.85 0.83],'GridAlpha',1);
exportgraphics(f,[fig 'F07_centreline.png'],'Resolution',140); close(f);

%% F08: profiles at stations
stx = [0 0.75 1.0 1.5 2.0 2.5]; f = figure('Position',[0 0 1500 560]); tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
qq = {S.U/U0, sqrt(S.uu)/U0, sqrt(S.vv)/U0, S.uv/U0^2}; qn = {'\langle u\rangle / U_0','u_{rms} / U_0','v_{rms} / U_0','\langle u''v''\rangle / U_0^2'};
yv = S.Y(:,1); hp = gobjects(numel(stx),1);
for q=1:4
  nexttile; hold on; grid on; set(gca,'GridColor',[0.88 0.88 0.86],'GridAlpha',1);
  for i=1:numel(stx), [~,ix]=min(abs(xv-stx(i))); hp(i)=plot(qq{q}(:,ix),yv,'-','Color',C{i},'LineWidth',1.8); end
  yline(0.5,':','Color',ink2,'LineWidth',0.8); yline(-0.5,':','Color',ink2,'LineWidth',0.8);
  xlabel(qn{q}); if q==1, ylabel('y/D'); end; ylim(yl);
end
lg = legend(hp, arrayfun(@(s) sprintf('x/D = %.2f',xv(find(abs(xv-s)==min(abs(xv-s)),1))),stx,'UniformOutput',false), 'Box','off'); lg.Layout.Tile='east';
exportgraphics(f,[fig 'F08_profiles.png'],'Resolution',140); close(f);

%% F09: instantaneous snapshots
f = figure('Position',[0 0 1500 1050]); tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for q=1:4
  us = fl(A.snap_u(:,:,q))/U0; vs = fl(A.snap_v(:,:,q))/U0; holes = isnan(us);
  m = ok; uf = us; vfl = vs; uf(~m)=0; vfl(~m)=0;
  uf = regionfill(uf, holes & m); vfl = regionfill(vfl, holes & m);
  uf(~m)=NaN; vfl(~m)=NaN;
  [dvx,~] = gradient(vfl, X(1,2)-X(1,1), Y(2,1)-Y(1,1)); [~,duy] = gradient(uf, X(1,2)-X(1,1), Y(2,1)-Y(1,1));
  w = dvx-duy; w(~m)=NaN; uf(~m)=NaN; vfl(~m)=NaN;
  nexttile; set(gca,'Color',nodata); hold on; pcolor(X,Y,w); shading interp; colormap(gca,div); caxis([-12 12]);
  cb=colorbar; cb.Label.String='\omega_z D / U_0'; s=2;
  quiver(X(1:s:end,1:s:end),Y(1:s:end,1:s:end),uf(1:s:end,1:s:end)-0*1,vfl(1:s:end,1:s:end),1.2,'Color',[ink 0.75],'LineWidth',0.5);
  cyl(); axis image; xlim(xl); ylim(yl); xlabel('x/D'); ylabel('y/D');
  title(sprintf('Snapshot #%d  (POD phase %.0f°)', A.snaps(q), A.phi(A.snaps(q))*180/pi));
end
exportgraphics(f,[fig 'F09_instantaneous.png'],'Resolution',140); close(f);

%% F10: POD energy and phase portrait
f = figure('Position',[0 0 1500 500]); tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; en = A.en_all; b = bar(1:20,100*en(1:20),0.7,'FaceColor',C{1},'EdgeColor','none'); hold on
ylabel('energy fraction [%]'); xlabel('mode index'); title('POD modal energy'); box off
nexttile; plot(1:200,100*cumsum(en(1:200)),'-','Color',C{1}); xlabel('number of modes'); ylabel('cumulative energy [%]'); grid on; title('Cumulative energy'); ylim([0 100]);
set(gca,'GridColor',[0.88 0.88 0.86],'GridAlpha',1);
nexttile; a = A.a; l = A.lam50; hold on
scatter(a(~A.use,1)/sqrt(2*l(1)),a(~A.use,2)/sqrt(2*l(2)),8,[0.75 0.75 0.73],'filled');
scatter(a(A.use,1)/sqrt(2*l(1)),a(A.use,2)/sqrt(2*l(2)),8,C{1},'filled','MarkerFaceAlpha',0.45);
plot(cos(th),sin(th),'--','Color',ink2,'LineWidth',1); axis equal; box off
xlabel('a_1 / \surd(2\lambda_1)'); ylabel('a_2 / \surd(2\lambda_2)'); title('Phase portrait of modes 1-2 (gray: excluded, r < 0.5)');
exportgraphics(f,[fig 'F10_pod_energy.png'],'Resolution',140); close(f);

%% F11: POD modes 1-6 (vorticity of mode + vectors)
f = figure('Position',[0 0 1500 980]); tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for k=1:6
  mu_ = fl(A.modes(:,:,1,k)); mv_ = fl(A.modes(:,:,2,k)); mm = ~isnan(mu_);
  uf=mu_; vfl=mv_;
  [dvx,~]=gradient(vfl,X(1,2)-X(1,1),Y(2,1)-Y(1,1)); [~,duy]=gradient(uf,X(1,2)-X(1,1),Y(2,1)-Y(1,1)); w=dvx-duy; w(~mm)=NaN;
  w = w/max(abs(w(:)),[],'omitnan');
  nexttile; set(gca,'Color',nodata); hold on; pcolor(X,Y,w); shading interp; colormap(gca,div); caxis([-1 1]);
  s=3; quiver(X(1:s:end,1:s:end),Y(1:s:end,1:s:end),mu_(1:s:end,1:s:end),mv_(1:s:end,1:s:end),1.2,'Color',[ink 0.7],'LineWidth',0.5);
  cyl(); axis image; xlim(xl); ylim(yl); xlabel('x/D'); ylabel('y/D'); title(sprintf('POD mode %d  (%.1f %% energy)',k,100*en(k)));
end
cb = colorbar; cb.Layout.Tile='east'; cb.Label.String = 'normalised mode vorticity';
exportgraphics(f,[fig 'F11_pod_modes.png'],'Resolution',140); close(f);

%% F12: phase-averaged vorticity (12 phases)
% The sign of the POD phase is arbitrary. Tracked vortex cores move UPSTREAM as
% phi increases (s9: about -0.10 D per 30-deg bin), so physical time runs with
% DECREASING phi. Tiles and GIF frames are therefore ordered from large to small phi.
ordT = A.nbin:-1:1;
f = figure('Position',[0 0 1800 1060]); tl = tiledlayout(3,4,'TileSpacing','tight','Padding','compact');
title(tl,'Time runs left to right, top to bottom (\phi decreasing; the sign of the POD phase is arbitrary)','FontSize',11);
for b=ordT
  w = fl(A.Wph(:,:,b))*D/U0; up=fl(A.Uph(:,:,b))/U0; vp=fl(A.Vph(:,:,b))/U0;
  nexttile; set(gca,'Color',nodata); hold on; pcolor(X,Y,w); shading interp; colormap(gca,div); caxis([-8 8]);
  s=3; quiver(X(1:s:end,1:s:end),Y(1:s:end,1:s:end),up(1:s:end,1:s:end)-1*0,vp(1:s:end,1:s:end),1.2,'Color',[ink 0.65],'LineWidth',0.45);
  cyl(); axis image; xlim(xl); ylim(yl); set(gca,'FontSize',8);
  title(sprintf('\\phi = %d° (n = %d)', round(-180+(b-0.5)*360/A.nbin), A.nph(b)),'FontSize',9);
end
cb = colorbar; cb.Layout.Tile='east'; cb.Label.String='\langle\omega_z\rangle_\phi D / U_0';
exportgraphics(f,[fig 'F12_phase_averaged_vorticity.png'],'Resolution',130); close(f);

%% F12b: GIF animation of phase average
gif = [fig 'F12b_phase_animation.gif'];
f = figure('Position',[0 0 900 640],'Color','w');
for rep=1:1
for b=ordT
  clf; w = fl(A.Wph(:,:,b))*D/U0; up=fl(A.Uph(:,:,b))/U0; vp=fl(A.Vph(:,:,b))/U0;
  ax=axes; set(ax,'Color',nodata); hold on; pcolor(X,Y,w); shading interp; colormap(div); caxis([-8 8]); cb=colorbar; cb.Label.String='\omega_z D/U_0';
  up2=up; vp2=vp; up2(isnan(up2))=0; vp2(isnan(vp2))=0; hs=streamslice(X,Y,up2-0,vp2,1.6); set(hs,'Color',[ink 0.6],'LineWidth',0.6);
  cyl(); axis image; xlim(xl); ylim(yl); xlabel('x/D'); ylabel('y/D'); title(sprintf('Phase-averaged vorticity, \\phi = %d°  (time: \\phi decreasing)', round(-180+(b-0.5)*360/A.nbin)));
  drawnow; fr = getframe(f); [im,map] = rgb2ind(fr.cdata,256,'nodither');
  if b==ordT(1), imwrite(im,map,gif,'gif','LoopCount',inf,'DelayTime',0.25); else, imwrite(im,map,gif,'gif','WriteMode','append','DelayTime',0.25); end
end
end
close(f);

%% F13: spectra and autocorrelation
f = figure('Position',[0 0 1500 520]); tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
sp = A.spec; nexttile; hold on
for i=1:4, semilogy(sp.f,sp.Pv(:,i)/U0^2,'-','Color',C{i},'LineWidth',1.5); end
set(gca,'YScale','log'); xlabel('f / f_s  [cycles per sample]'); ylabel('PSD of v / U_0^2'); grid on; set(gca,'GridColor',[0.88 0.88 0.86],'GridAlpha',1);
legend(arrayfun(@(i) sprintf('(%.2f, %.2f)',sp.xy(i,1),sp.xy(i,2)),1:4,'UniformOutput',false),'Box','off','Location','southwest');
title('Spectra of v on the centreline (x/D, y/D)'); xlim([0 0.5]);
nexttile; hold on; semilogy(sp.f,sp.Pa1,'-','Color',C{1},'LineWidth',1.5); semilogy(sp.f,sp.Pa3,'-','Color',C{2},'LineWidth',1.5); set(gca,'YScale','log');
xlabel('f / f_s'); ylabel('PSD'); legend({'a_1','a_3'},'Box','off','Location','southwest'); title('Spectra of POD coefficients'); grid on; xlim([0 0.5]);
set(gca,'GridColor',[0.88 0.88 0.86],'GridAlpha',1);
nexttile; hold on; stem(1:8,A.ac,'filled','Color',C{1},'MarkerSize',5,'LineWidth',1.5); yline(0,'-','Color',ink2,'LineWidth',0.8);
xlabel('lag [samples]'); ylabel('autocorrelation of a_1'); title('Sign flips every sample: shedding near f_s/2'); ylim([-1 1]); box off
exportgraphics(f,[fig 'F13_spectra.png'],'Resolution',140); close(f);

%% F14: PDFs and convergence
cv = A.conv; f = figure('Position',[0 0 1500 500]); tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; hold on; n = (1:numel(cv.mean))';
plot(n,cv.mean/U0,'-','Color',C{1}); plot(n,cv.rms/U0,'-','Color',C{2});
xlabel('number of samples'); ylabel('running statistic / U_0'); legend({'mean u','u_{rms}'},'Box','off','Location','east');
title(sprintf('Convergence at (x/D, y/D) = (%.2f, %.2f)',cv.pt)); grid on; set(gca,'GridColor',[0.88 0.88 0.86],'GridAlpha',1);
nexttile; hold on; histogram(cv.pdf_u/U0,60,'Normalization','pdf','FaceColor',C{1},'EdgeColor','w','FaceAlpha',0.85);
histogram(cv.pdf_v/U0,60,'Normalization','pdf','FaceColor',C{2},'EdgeColor','w','FaceAlpha',0.6);
xlabel('velocity / U_0'); ylabel('PDF'); legend({'u','v'},'Box','off'); title(sprintf('PDF in the shear layer (%.2f, %.2f)',cv.pt));
nexttile; hold on; g1 = cv.pdf_u_cl(~isnan(cv.pdf_u_cl)); g2 = cv.pdf_v_cl(~isnan(cv.pdf_v_cl));
histogram(g1/U0,60,'Normalization','pdf','FaceColor',C{1},'EdgeColor','w','FaceAlpha',0.85);
histogram(g2/U0,60,'Normalization','pdf','FaceColor',C{2},'EdgeColor','w','FaceAlpha',0.6);
xlabel('velocity / U_0'); legend({'u','v'},'Box','off'); title(sprintf('PDF on the centreline (%.2f, %.2f)',cv.pt_cl));
exportgraphics(f,[fig 'F14_pdf_convergence.png'],'Resolution',140); close(f);

%% F15: bulk (empty flume) profiles + comparison
bp = A.bprof; f = figure('Position',[0 0 1500 560]); tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; hold on; grid on; set(gca,'GridColor',[0.88 0.88 0.86],'GridAlpha',1);
fill([bp.U-bp.Ustd_x; flipud(bp.U+bp.Ustd_x)],[bp.y; flipud(bp.y)],C{1},'FaceAlpha',0.15,'EdgeColor','none');
p1=plot(bp.U,bp.y,'-','Color',C{1});
[~,ixu] = min(abs(xv-(-0.85))); p2=plot(S.U(:,ixu),yv,'-','Color',C{2});
[~,ixd] = min(abs(xv-2.5)); p3=plot(S.U(:,ixd),yv,'-','Color',C{3});
xlabel('\langle u\rangle [m/s]'); ylabel('y/D'); yline([-0.5 0.5],':','Color',ink2,'LineWidth',0.8,'HandleVisibility','off');
legend([p1 p2 p3],{'empty flume (x-averaged)',sprintf('with cylinder, x/D = %.2f',xv(ixu)),sprintf('with cylinder, x/D = %.2f',xv(ixd))},'Box','off','Location','northwest');
title('Approach flow vs. flow with cylinder');
nexttile; hold on; grid on; set(gca,'GridColor',[0.88 0.88 0.86],'GridAlpha',1);
plot(bp.urms,bp.y,'-','Color',C{1}); plot(bp.vrms,bp.y,'-','Color',C{2}); xlabel('rms [m/s]'); legend({'u_{rms}','v_{rms}'},'Box','off'); title('Empty flume fluctuations');
nexttile; set(gca,'Color',nodata); hold on; F = fl(A.Ub); pcolor(X,Y,F); shading interp; colormap(gca,seq); cb=colorbar; cb.Label.String='\langle u\rangle [m/s]';
axis image; xlim(xl); ylim(yl); xlabel('x/D'); title(sprintf('Empty flume mean field (N = %d)',A.nb_t));
plot(0.5*cos(th),0.5*sin(th),'--','Color',C{2},'LineWidth',1);
exportgraphics(f,[fig 'F15_bulk_flow.png'],'Resolution',140); close(f);

%% F16: calibration
cal = double(imread([root 'Calibration/Generated000000.T000.D000.P000.H000.LA.TIF']));
f = figure('Position',[0 0 1500 600]); tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; imagesc(cal); colormap(gca,gray); caxis(prctile(cal(:),[1 99.7])); axis image; hold on
plot(386.9+188.45*cos(th),588.5+188.45*sin(th),'-','Color',C{2},'LineWidth',1.2); plot(386.9,588.5,'+','Color',C{2},'MarkerSize',14);
title('Calibration image: circle fit to the illuminated cylinder rim'); xlabel('x [px]'); ylabel('row [px]');
nexttile; imagesc(cal); colormap(gca,gray); caxis(prctile(cal(:),[1 99.9])); axis image; hold on; xlim([150 330]); ylim([420 760]);
plot(386.9+186.0*cos(th),588.5+186.0*sin(th),':','Color',C{1},'LineWidth',1.2);
plot(386.9+188.45*cos(th),588.5+188.45*sin(th),'-','Color',C{2},'LineWidth',1.2);
plot(386.9+191.3*cos(th),588.5+191.3*sin(th),':','Color',C{3},'LineWidth',1.2);
legend({'inner half-max D=372 px','rim peak D=376.9 px','outer half-max D=382.6 px'},'Box','off','TextColor','w','Location','southeast');
title('Zoom: rim edge definitions (uncertainty band)');
exportgraphics(f,[fig 'F16_calibration.png'],'Resolution',140); close(f);
disp('figures done');

function h = tiledlayout_handle_fallback(f)
  h = findobj(f,'Type','tiledlayout'); h = h(1);
end
