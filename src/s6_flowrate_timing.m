% Step 6: discharge Q from the empty-flume profile + check of an assumed 1 s frame interval
clear; res = piv_path('results');
fig = [res 'figures/'];
A = load([res 'analysis.mat']); S = load([res 'cyl_stats.mat'],'D','nu','rho','Utop');
D = S.D; nu = S.nu; rho = S.rho; b = 0.060; h = 0.118;   % width 60 mm (user, 2026-11-21); worksheet entry 0.60 m was wrong by 10x
%% ---------- 1. measured part of the profile ----------
y = A.bprof.y*D; u = A.bprof.U; g = ~isnan(u); y = y(g); u = u(g); [y,o] = sort(y); u = u(o);
ybot = min(y); ytop = max(y); Lfov = ytop - ybot;
Ucols = A.Ub(:,3:85); Ucols = Ucols(g,:); Ucols = Ucols(o,:); Ucols = fillmissing(Ucols,'linear',1,'EndValues','nearest'); Ucols = Ucols(:,all(~isnan(Ucols),1)); ncol = size(Ucols,2);
I_fov = trapz(y,u);                                     % m^2/s, measured part
u_bot = u(1); u_top = u(end);
pt = polyfit(y(end-5:end), u(end-5:end), 1); slope_top = pt(1);     % trend of the top 6 rows
fprintf('FOV covers %.1f mm of %.0f mm depth (%.0f %%), measured integral %.5f m^2/s, mean %.4f m/s\n', Lfov*1e3, h*1e3, 100*Lfov/h, I_fov, I_fov/Lfov);
%% ---------- 2. Monte Carlo over the unknowns ----------
rng(1); N = 2e5;
yb_lo_feas = ytop - h; yb_hi_feas = ybot;                % bed must lie below the FOV, surface above it
% (a) conservative: bed anywhere in the feasible range
% (b) best estimate: velocity maximum at y = 39 mm is read as the velocity dip at z_max/h = 0.80-0.90
[~,im] = max(u); ymax = y(im);
yb_lo_dip = ymax - 0.90*h; yb_hi_dip = ymax - 0.80*h;
% (c) given by the user: cylinder centre about 60 mm above the bed -> y_b = -60 mm, reading uncertainty +-3 mm
cases = {'conservative', yb_lo_feas, yb_hi_feas; 'best', -0.063, -0.057; 'dip', max(yb_lo_dip,yb_lo_feas), min(yb_hi_dip,yb_hi_feas)};
Q = struct();
for c = 1:3
  yb = cases{c,2} + (cases{c,3}-cases{c,2})*rand(N,1);
  gb = ybot - yb;  gt = h - Lfov - gb;                         % bottom / top gaps
  mb = 3 + 7*rand(N,1);                                       % near-bed power law 1/3 .. 1/10 (widened after review)
  Ibot = u_bot*gb.*mb./(mb+1);
  wt = rand(N,1);                                             % top: between constant and linear trend
  Itop = u_top*gt + wt.*0.5*slope_top.*gt.^2;
  colpick = randi(ncol,N,1);                                  % statistical: resample columns (bootstrap of x-positions)
  Ifov_bs = arrayfun(@(k) trapz(y,Ucols(:,k)), 1:ncol)'; Ifov_bs = Ifov_bs - mean(Ifov_bs) + I_fov;
  Ifv = I_fov + (Ifov_bs(colpick)-I_fov)/sqrt(ncol);
  calib = 1 + 0.015*(2*rand(N,1)-1);                          % magnification +-1.5 %
  plock = 1 + 0.01*(2*rand(N,1)-1);                           % peak-locking bias of the mean +-1 %
  dl = min((0.5 + 0.5*rand(N,1))*h, b/2); mw = 5 + 4*rand(N,1);   % side-wall layers; in a 60 mm channel they meet at the centre plane (dl = b/2 -> kw = m/(m+1))
  kw = 1 - 2*dl./(b*(mw+1));
  ub_c = (Ibot + Ifv + Itop)/h .* calib .* plock;             % centre-plane depth-mean velocity
  Qc = b*h*ub_c.*kw;
  Q.(cases{c,1}) = struct('Q',Qc,'ub',Qc/(b*h),'ubc',ub_c,'kw',kw,'yb',yb);
  pr = prctile(Qc,[2.5 50 97.5])*1e3;
  fprintf('%-12s: y_b in [%.1f, %.1f] mm -> Q = %.2f l/s (95%%: %.2f .. %.2f), u_b = %.4f m/s, +-%.1f %%\n', cases{c,1}, cases{c,2}*1e3, cases{c,3}*1e3, ...
     mean(Qc)*1e3, pr(1), pr(3), mean(Qc)/(b*h), 100*(pr(3)-pr(1))/2/pr(2));
end
% variance budget for the best case (one factor varied at a time, others at mid values)
Qb = Q.best; ub0 = mean(Qb.ub);
fprintf('naive Q = b*y0*mean(FOV) = %.2f l/s\n', b*h*I_fov/Lfov*1e3);
% sensitivity: one-at-a-time ranges (best case)
ybm = mean([cases{2,2:3}]); gbm = ybot-ybm; gtm = h-Lfov-gbm;
f = @(gb,mb,wt,cal,pl,dl,mw) b*(u_bot*gb*mb/(mb+1) + I_fov + u_top*(h-Lfov-gb) + wt*0.5*slope_top*(h-Lfov-gb)^2)*cal*pl*(1-2*dl/(b*(mw+1)));
q0 = f(gbm,7,0.5,1,1,min(0.75*h,b/2),7);
sens = { 'bed position', f(ybot-cases{2,3},7,0.5,1,1,min(0.75*h,b/2),7), f(ybot-cases{2,2},7,0.5,1,1,min(0.75*h,b/2),7);
         'near-bed law', f(gbm,4,0.5,1,1,min(0.75*h,b/2),7), f(gbm,10,0.5,1,1,min(0.75*h,b/2),7);
         'top extrapolation', f(gbm,7,0,1,1,min(0.75*h,b/2),7), f(gbm,7,1,1,1,min(0.75*h,b/2),7);
         'magnification', f(gbm,7,0.5,0.985,1,min(0.75*h,b/2),7), f(gbm,7,0.5,1.015,1,min(0.75*h,b/2),7);
         'peak locking', f(gbm,7,0.5,1,0.99,min(0.75*h,b/2),7), f(gbm,7,0.5,1,1.01,min(0.75*h,b/2),7);
         'side-wall layer', f(gbm,7,0.5,1,1,min(h,b/2),5), f(gbm,7,0.5,1,1,min(0.5*h,b/2),9)};
fprintf('central Q = %.3f l/s; one-at-a-time half ranges:\n', q0*1e3);
for k=1:size(sens,1), fprintf('  %-18s %.3f .. %.3f l/s  (+-%.1f %%)\n', sens{k,1}, min(sens{k,2:3})*1e3, max(sens{k,2:3})*1e3, 50*abs(sens{k,3}-sens{k,2})/q0); end
%% ---------- 3. derived numbers with u_b ----------
ub = mean(Q.best.ub); Qbest = mean(Q.best.Q);
mu = rho*nu; dp = 50e-6; rp = 1016; gg = 9.81;
tau_p = rp*dp^2/(18*mu); ws = (rp-rho)*gg*dp^2/(18*mu);
eps_ = ub^3/D; tau_eta = sqrt(nu/eps_); eta_K = (nu^3/eps_)^0.25; v_K = (nu*eps_)^0.25;
nd2 = struct('Q',Qbest,'Q_lo',prctile(Q.best.Q,2.5),'Q_hi',prctile(Q.best.Q,97.5),'Q_cons_lo',prctile(Q.conservative.Q,2.5),'Q_cons_hi',prctile(Q.conservative.Q,97.5), ...
   'ub',ub,'ReD',ub*D/nu,'Re_h',ub*h/nu,'Fr',ub/sqrt(gg*h),'beta',D/h,'tau_p',tau_p,'ws',ws,'eps',eps_,'tau_eta',tau_eta,'eta_K',eta_K,'v_K',v_K,'St_k',tau_p/tau_eta);
disp(nd2)
%% ---------- 4. assumed frame interval 1 s ----------
fs = 1.0;
a = A.a; lam = A.lam50; n = (0:size(a,1)-1)';
z = (a(:,1)/sqrt(2*lam(1)) + 1i*a(:,2)/sqrt(2*lam(2))) .* (-1).^n;   % demodulate by f_s/2
% heavily averaged Welch spectrum (short segments) -> smooth hump; use its centroid and spread
nfft = 128; w = 0.5-0.5*cos(2*pi*(0:nfft-1)'/(nfft-1)); nseg = floor((numel(z)-nfft)/(nfft/2))+1;
nz = 32*nfft; P = zeros(nz,1);
for k=1:nseg, seg = (z((k-1)*nfft/2+(1:nfft)) - mean(z)).*w; F = fft(seg,nz); P = P + abs(F).^2; end
P = fftshift(P)/nseg; df = ((-nz/2):(nz/2-1))'/nz;               % offset delta from f_s/2, cycles/sample
win = abs(df) <= 0.1; Pw = P(win) - min(P(win)); dw = df(win);
delta = sum(dw.*Pw)/sum(Pw);                                    % centroid
bw = 2*sqrt(sum((dw-delta).^2.*Pw)/sum(Pw));                    % 2 x rms spread (includes window broadening)
% phase-step estimate as an independent check
dphi = angle(exp(1i*(diff(A.phi)))); dphi = dphi(A.use(1:end-1) & A.use(2:end));
delta_ph = abs(angle(mean(exp(1i*dphi))))/(2*pi) - 0.5;
fprintf('phase-step estimate: |delta| = %.4f cycles/sample\n', abs(delta_ph));
fA = (0.5 + abs(delta))*fs; fB = (0.5 - abs(delta))*fs;
fprintf('demodulated spectrum centroid |delta| = %.4f cycles/sample, 2-sigma spread %.4f\n', abs(delta), bw);
fprintf('f_s = 1 Hz -> f_shed = %.3f or %.3f Hz (alias ambiguity), period %.2f s\n', fA, fB, 1/fA);
fsh = [0.4781 0.4781]*fs;   % phase-step estimate with the sign fixed by vortex travel (s9); was [fB fA] from the centroid
St_ub = fsh*D/ub; St_U0 = fsh*D/A.U0; St_gap = fsh*D/S.Utop; St_blk = St_ub*(1-D/h);
fprintf('St based on u_b: %.3f .. %.3f; on U0: %.3f .. %.3f; on gap velocity: %.3f .. %.3f; blockage-corrected (u_b/(1-beta)): %.3f .. %.3f\n', St_ub, St_U0, St_gap, St_blk);
fprintf('St spread from spectral width: +-%.3f\n', 0.5*bw*fs*D/ub);
% other aliases for f_s = 1 Hz
for k = 1:2, fprintf('alias k=%d: f = %.2f Hz -> St(u_b) = %.2f\n', k, (k+0.5)*fs, (k+0.5)*fs*D/ub); end
% which f_s would be needed for St(u_b) = 0.21 (k = 0,1,2)?
for k=0:2, fprintf('St=0.21 would need f_s = %.3f Hz (k=%d)\n', 0.21*ub/D/(k+0.5), k); end
fprintf('record length at 1 Hz: cyl %.1f min, bulk %.1f min; shedding cycles in record: %.0f\n', size(a,1)/60, A.nb_t/60, size(a,1)*0.5);
T = struct('fs',fs,'delta',delta,'bw',bw,'fsh',fsh,'St_ub',St_ub,'St_U0',St_U0,'St_gap',St_gap,'St_blk',St_blk,'P',P,'df',df);
save([res 'flowrate_timing.mat'],'Q','nd2','T','y','u','ybot','ytop','cases','sens','q0','u_bot','u_top','slope_top','I_fov','Lfov','-v7.3');

%% ---------- 5. figures ----------
hex = @(hh) sscanf(hh(2:end),'%2x%2x%2x',[1 3])/255;
C = cellfun(hex,{'#2a78d6','#eb6834','#1baf7a','#eda100','#e87ba4'},'UniformOutput',false); ink2 = hex('#52514e');
set(groot,'defaultAxesFontSize',10,'defaultAxesTickDir','out','defaultAxesBox','off','defaultFigureColor','w','defaultLineLineWidth',2,'defaultAxesTitleFontWeight','normal');
% F17: discharge construction
fg = figure('Position',[0 0 1500 560]); tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; hold on; grid on; set(gca,'GridColor',[0.88 0.88 0.86],'GridAlpha',1);
ybm = mean([cases{2,2:3}]); gb = ybot - ybm; zz = linspace(0,1,40)';
patch([0 0 0.09 0.09]*1, [cases{1,2} cases{1,3} cases{1,3} cases{1,2}]*1e3, [0.85 0.85 0.83],'EdgeColor','none','FaceAlpha',0.6);
patch([0 0 0.09 0.09], [cases{2,2} cases{2,3} cases{2,3} cases{2,2}]*1e3, C{2},'EdgeColor','none','FaceAlpha',0.25);
p1 = plot(u, y*1e3, '-', 'Color', C{1});
p2 = plot(u_bot*zz.^(1/7), (ybm + zz*gb)*1e3, '--', 'Color', C{1});
gt = h - Lfov - gb; yt = linspace(ytop, ytop+gt, 20)'; p3 = plot(u_top + 0.5*slope_top*(yt-ytop)*0 + 0*yt, yt*1e3, ':', 'Color', C{1});
yline(ybm*1e3,'-','Color',C{2},'LineWidth',1.2); yline((ybm+h)*1e3,'-','Color',C{3},'LineWidth',1.2);
text(0.002,ybm*1e3+2,'bed (given: 60 mm below centre)','Color',ink2,'FontSize',9); text(0.002,(ybm+h)*1e3-3,'free surface','Color',ink2,'FontSize',9);
text(0.002,cases{1,3}*1e3-4,'gray: bed geometrically possible  |  orange: 60 ± 3 mm','Color',ink2,'FontSize',8.5);
xlim([0 0.09]); ylim([cases{1,2}*1e3-5, (ybm+h)*1e3+6]); xlabel('\langle u\rangle [m/s]'); ylabel('y relative to cylinder centre [mm]');
legend([p1 p2 p3],{'measured (empty flume)','near-bed power law','top extrapolation'},'Box','off','Location','west');
title('Depth profile used for the discharge');
nexttile; hold on;
histogram(Q.conservative.Q*1e3, 80,'Normalization','pdf','FaceColor',[0.7 0.7 0.68],'EdgeColor','none');
histogram(Q.best.Q*1e3, 80,'Normalization','pdf','FaceColor',C{2},'EdgeColor','none','FaceAlpha',0.8);
xline(b*h*I_fov/Lfov*1e3,'--','Color',ink2,'LineWidth',1.2);
legend({'bed anywhere feasible','bed 60 ± 3 mm below centre','naive b\cdoty_0\cdotmean(FOV)'},'Box','off','Location','northwest');
xlabel('Q [l/s]'); ylabel('PDF'); title('Monte Carlo distribution of Q'); box off
nexttile; hold on; names = sens(:,1); lo = (cell2mat(sens(:,2))-q0)/q0*100; hi = (cell2mat(sens(:,3))-q0)/q0*100;
[~,ord] = sort(abs(hi-lo));
for k=1:numel(ord), i = ord(k); plot([min(lo(i),hi(i)) max(lo(i),hi(i))],[k k],'-','Color',C{1},'LineWidth',8); end
xline(0,'-','Color',ink2,'LineWidth',0.8); set(gca,'YTick',1:numel(ord),'YTickLabel',names(ord)); xlabel('change of Q [%]');
title(sprintf('Sensitivity (central Q = %.2f l/s)', q0*1e3)); grid on; set(gca,'GridColor',[0.88 0.88 0.86],'GridAlpha',1); ylim([0.4 numel(ord)+0.6]);
exportgraphics(fg,[fig 'F17_discharge.png'],'Resolution',140); close(fg);
% F18: timing check
fg = figure('Position',[0 0 1500 520]); tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; sp = A.spec; hold on
semilogy(sp.f*fs, sp.Pa1, '-', 'Color', C{1}, 'LineWidth', 1.5); set(gca,'YScale','log');
xline(fA,'--','Color',C{2},'LineWidth',1.2); xlabel('f [Hz]  (f_s = 1 Hz assumed)'); ylabel('PSD of a_1'); xlim([0 0.5]);
title('a_1 spectrum on an assumed 1 Hz time axis'); grid on; set(gca,'GridColor',[0.88 0.88 0.86],'GridAlpha',1);
nexttile; plot(df*fs*1e3, P/max(P), '-', 'Color', C{1}); hold on; xline(delta*fs*1e3,'--','Color',C{2},'LineWidth',1.2);
xlim([-100 100]); xlabel('offset from f_s/2 [mHz]'); ylabel('normalised PSD'); title(sprintf('Demodulated a_1 + i a_2: centroid %+.1f mHz, spread ±%.1f mHz', delta*1e3, bw*1e3));
grid on; set(gca,'GridColor',[0.88 0.88 0.86],'GridAlpha',1);
nexttile; hold on; cats = {'u_b','u_b/(1-D/y_0)','U_0 (bulk at cyl.)','gap velocity'}; vals = [mean(St_ub) mean(St_blk) mean(St_U0) mean(St_gap)];
patch([0.18 0.24 0.24 0.18],[0.4 0.4 4.6 4.6],C{3},'EdgeColor','none','FaceAlpha',0.15);
for k=1:4, plot(vals(k),k,'o','MarkerFaceColor',C{1},'MarkerEdgeColor','w','MarkerSize',9); text(vals(k)+0.01,k,sprintf('%.3f',vals(k)),'FontSize',9,'Color',ink2); end
set(gca,'YTick',1:4,'YTickLabel',cats); xlim([0.1 0.35]); ylim([0.4 4.6]); xlabel('Strouhal number = f D / U_{ref}');
title('Strouhal number for f_s = 1 Hz (band: 0.18-0.24)'); grid on; set(gca,'GridColor',[0.88 0.88 0.86],'GridAlpha',1);
exportgraphics(fg,[fig 'F18_timing_check.png'],'Resolution',140); close(fg);
disp('done');
