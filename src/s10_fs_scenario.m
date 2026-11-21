% Step 10: what changes if the frame rate is f_s = 2 Hz instead of the assumed 1 Hz?
% Only quantities that carry seconds depend on f_s. Everything built from single
% fields (means, rms, profiles, Lr, POD modes, phase averages, Gamma, epsilon, Q)
% uses only the pulse separation dt = 3000 us and the magnification, so it is unchanged.
clear; res = piv_path('results'); fig = [res 'figures/'];
A = load([res 'analysis.mat'],'U0'); S = load([res 'cyl_stats.mat'],'D','nu','Utop');
L = load([res 'literature.mat'],'lambda_mean','lambda_pod');
D = S.D; U0 = A.U0; h = 0.118;
ub = 0.0542;                           % b = 60 mm, narrow-channel side-wall correction (s6)
Ublk = ub/(1 - D/h); Utop = S.Utop;
frac = 0.4781;                         % f/f_s = n + frac, sign fixed by vortex travel (s9)
lam_core = L.lambda_mean; lam_pod = L.lambda_pod;
Ncyl = 4069; Tbulk_1Hz = 11.5*60;
refs = {'u_b', ub; 'U_0', U0; 'blockage-corrected', Ublk; 'gap (above cyl.)', Utop};
fprintf('ub %.4f  U0 %.4f  Ublk %.4f  Utop %.4f m/s; lambda core %.2f D, POD %.2f D\n', ub, U0, Ublk, Utop, lam_core, lam_pod);
for fs = [1 2]
  f = frac*fs; fprintf('\n===== f_s = %g Hz, n = 0 =====\n', fs);
  fprintf('f_shed = %.3f Hz, period %.2f s (%.2f samples)\n', f, 1/f, 1/frac);
  for k = 1:4, fprintf('St (%s) = %.3f\n', refs{k,1}, f*D/refs{k,2}); end
  fprintf('n = 1: f = %.2f Hz, St_blk = %.2f\n', (1+frac)*fs, (1+frac)*fs*D/Ublk);
  fprintf('Uc/U0 vortex cores = %.2f, phase velocity (POD) = %.2f\n', lam_core*D*f/U0, lam_pod*D*f/U0);
  fprintf('record: cyl %.1f min, bulk %.2f min, cycles %.0f\n', Ncyl/fs/60, Tbulk_1Hz/fs/60, Ncyl*frac);
  fprintf('amplitude e-fold %.2f s (0.88 cycles), a1 envelope 1/e %.1f s (3 cycles)\n', 1.85/fs, 6/fs);
end
% velocity scale that 2 Hz would need for St_blk = 0.21
fprintf('\nfor St_blk = 0.21 at 2 Hz the blockage-corrected velocity would have to be %.3f m/s (measured %.3f, x%.2f)\n', ...
  2*frac*D/0.21, Ublk, 2*frac*D/0.21/Ublk);
fprintf('f_s window for St_blk in [0.18 0.24], n = 0: %.2f - %.2f Hz\n', 0.18*Ublk/(frac*D), 0.24*Ublk/(frac*D));

%% figure: sensitivity to f_s
hex = @(s) sscanf(s(2:end),'%2x%2x%2x')'/255;
C = cellfun(hex,{'#2a78d6','#eb6834','#1baf7a','#eda100'},'UniformOutput',false);
ink = hex('#52514e'); grid_c = [0.88 0.88 0.86];
fsv = linspace(0.3,3,400);
f1 = figure('Units','centimeters','Position',[2 2 24 9.5],'Color','w');
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on
patch([fsv(1) fsv(end) fsv(end) fsv(1)],[0.18 0.18 0.24 0.24],hex('#1baf7a'),'FaceAlpha',0.13,'EdgeColor','none','DisplayName','plausible Strouhal number 0.18-0.24');
yline(0.215,':','Color',ink,'LineWidth',1,'DisplayName','Tremblay exp. 0.215');
plot(fsv, frac*fsv*D/Ublk,'-','Color',C{1},'LineWidth',2,'DisplayName','n = 0');
plot(fsv, (1+frac)*fsv*D/Ublk,'--','Color',C{1},'LineWidth',1.5,'DisplayName','n = 1');
xw = [0.18 0.24]*Ublk/(frac*D);
xline(xw,'-','Color',hex('#1baf7a'),'LineWidth',0.8,'HandleVisibility','off');
for fs = [1 2]
  s = frac*fs*D/Ublk; c = C{1+(fs==2)};
  plot(fs, s,'o','MarkerSize',8,'MarkerFaceColor',c,'MarkerEdgeColor','w','LineWidth',1.5,'HandleVisibility','off');
  if fs == 1, text(fs+0.08, s-0.045, sprintf('%g Hz: Strouhal number = %.3f', fs, s),'Color',ink,'FontSize',9);
  else, text(fs+0.1, s-0.06, {sprintf('%g Hz:', fs), sprintf('Strouhal number = %.3f', s)},'Color',ink,'FontSize',9); end
end
xlim([0.3 3]); ylim([0 0.8]); xlabel('assumed frame rate f_s [Hz]'); ylabel('Strouhal number = f D / U_{blk}');
title(sprintf('(a) Strouhal number vs f_s; window %.2f-%.2f Hz', xw)); grid on; set(gca,'GridColor',grid_c,'GridAlpha',1);
legend('Box','off','Location','southoutside','NumColumns',2);
nexttile; hold on
yline(1,'-','Color',ink,'LineWidth',0.8,'DisplayName','U_c = U_0');
plot(fsv, lam_core*D*frac*fsv/U0,'-','Color',C{1},'LineWidth',2,'DisplayName',sprintf('vortex cores (\\lambda = %.2f D)',lam_core));
plot(fsv, lam_pod*D*frac*fsv/U0,'-','Color',C{3},'LineWidth',2,'DisplayName',sprintf('POD phase velocity (\\lambda = %.2f D)',lam_pod));
for fs = [1 2]
  for q = 1:2
    lam = [lam_core lam_pod]; v = lam(q)*D*frac*fs/U0; c = C{[1 3]}; cc = {C{1},C{3}};
    plot(fs, v,'o','MarkerSize',8,'MarkerFaceColor',cc{q},'MarkerEdgeColor','w','LineWidth',1.5,'HandleVisibility','off');
    text(fs+0.07, v, sprintf('%.2f',v),'Color',ink,'FontSize',9);
  end
end
xline([1 2],':','Color',ink,'HandleVisibility','off');
xlim([0.3 3]); ylim([0 2]); xlabel('assumed frame rate f_s [Hz]'); ylabel('U_c / U_0  (n = 0)');
title('(b) convection speed in the formation region'); grid on; set(gca,'GridColor',grid_c,'GridAlpha',1);
legend('Box','off','Location','northwest');
exportgraphics(f1,[fig 'N07_fs_sensitivity.png'],'Resolution',140); close(f1);
