function L = lit_tremblay()
% Reference data from F. Tremblay (2001), "Direct and large-eddy simulation of flow
% around a circular cylinder at subcritical Reynolds numbers", Diss. TU Muenchen.
% Re = 3900, U0 = free-stream velocity, x measured from the cylinder centre.
% Curves DIGITISED BY EYE from the printed figures (accuracy about +-0.02 in the
% ordinate, +-0.05 D in position) -- adequate for trend comparison, not for fine detail.
% DNS = present DNS of the thesis; LS = PIV of Lourenco & Shih (exp.);
% OW = hot wire of Ong & Wallace (exp.); Ma = DNS of Ma et al. (Case I).

% ---- Table 3.1 (mean flow parameters) -- values printed in the thesis ----
L.tab.name = {'Exp.','Ma et al. DNS','Tremblay DNS'};
L.tab.CD   = [0.99 NaN 1.03];
L.tab.Cpb  = [-0.88 -0.96 -0.93];
L.tab.sep  = [86 NaN 85.7];
L.tab.Lr   = [1.4 1.12 1.30];      % from the REAR surface, in D
L.tab.dLr  = [0.1 0 0];
L.tab.St   = [0.215 0.203 0.220];
L.tab.dSt  = [0.005 0 0];
L.Re = 3900;

% ---- Fig. 3.15: centreline mean streamwise velocity ----
L.cl.DNS = [0.5 0; 0.7 -0.08; 0.9 -0.17; 1.1 -0.23; 1.25 -0.25; 1.4 -0.22; 1.6 -0.12; 1.8 0.0; ...
            2.0 0.14; 2.2 0.30; 2.4 0.44; 2.6 0.56; 2.8 0.64; 3.0 0.69; 3.5 0.745; 4.0 0.765; ...
            5 0.775; 6 0.78; 8 0.78; 10 0.78];
L.cl.LS  = [0.55 -0.02; 0.7 -0.08; 0.9 -0.18; 1.05 -0.23; 1.2 -0.25; 1.35 -0.23; 1.5 -0.16; ...
            1.6 -0.07; 1.7 0.03; 1.8 0.10; 1.9 0.18; 2.0 0.24; 2.2 0.30; 2.4 0.37; 2.6 0.42; ...
            2.8 0.44; 3.0 0.47; 3.2 0.55; 3.4 0.65; 3.6 0.70; 3.8 0.72; 4.0 0.74; 4.2 0.77];
L.cl.OW  = [4.0 0.75; 5.0 0.77; 6.0 0.79; 7.0 0.775; 10 0.78];

% ---- Figs. 3.16/3.17/4.38: mean u/U0 profiles (upper half, symmetric) ----
% columns: y/D, u/U0
L.pu.x = [1.06 1.54 2.02 3.00];
L.pu.DNS{1} = [0 -0.23; 0.2 -0.17; 0.35 -0.02; 0.45 0.20; 0.55 0.60; 0.62 1.00; 0.70 1.32; 1.0 1.27; 1.5 1.19; 2.0 1.12; 3.0 1.05];
L.pu.LS{1}  = [0 -0.22; 0.1 -0.20; 0.2 -0.15; 0.3 -0.05; 0.4 0.10; 0.45 0.25; 0.5 0.45; 0.6 0.75; 0.65 1.05; 0.7 1.20; 1.0 1.22; 1.5 1.18];
L.pu.DNS{2} = [0 -0.11; 0.2 0.02; 0.35 0.30; 0.45 0.52; 0.55 0.75; 0.7 1.00; 0.85 1.10; 1.0 1.12; 1.5 1.10; 2.0 1.08; 3.0 1.05];
L.pu.LS{2}  = [0 -0.13; 0.2 0.00; 0.35 0.27; 0.45 0.48; 0.55 0.70; 0.7 0.95; 0.85 1.07; 1.0 1.10; 1.5 1.10];
L.pu.DNS{3} = [0 0.30; 0.2 0.38; 0.4 0.58; 0.5 0.68; 0.6 0.80; 0.8 0.97; 1.0 1.04; 1.5 1.07; 2.0 1.06; 3.0 1.05];
L.pu.LS{3}  = [0 0.22; 0.2 0.30; 0.4 0.50; 0.5 0.62; 0.6 0.74; 0.8 0.92; 1.0 1.02; 1.5 1.07];
L.pu.DNS{4} = [0 0.68; 0.2 0.72; 0.4 0.80; 0.55 0.86; 0.8 0.95; 1.0 1.00; 1.5 1.04; 2.0 1.05; 3.0 1.05];
L.pu.LS{4}  = [0 0.52; 0.2 0.57; 0.4 0.70; 0.55 0.80; 0.8 0.92; 1.0 0.98; 1.5 1.07];

% ---- Figs. 3.22/3.23: <u'u'>/U0^2 profiles (upper half) ----
L.uu.DNS{1} = [0 0.10; 0.2 0.12; 0.3 0.16; 0.45 0.16; 0.55 0.20; 0.60 0.26; 0.65 0.12; 0.75 0.03; 1.0 0.005; 1.5 0];
L.uu.LS{1}  = [0 0.04; 0.2 0.08; 0.3 0.14; 0.4 0.17; 0.5 0.20; 0.55 0.23; 0.65 0.10; 0.75 0.03; 1.0 0.01; 1.5 0];
L.uu.DNS{2} = [0 0.10; 0.15 0.14; 0.3 0.24; 0.4 0.27; 0.5 0.22; 0.6 0.12; 0.8 0.03; 1.0 0.01; 1.5 0];
L.uu.LS{2}  = [0 0.08; 0.15 0.10; 0.3 0.16; 0.45 0.18; 0.55 0.15; 0.7 0.06; 0.9 0.02; 1.2 0.005; 1.5 0];
L.uu.DNS{3} = [0 0.09; 0.15 0.11; 0.3 0.14; 0.5 0.16; 0.7 0.12; 0.9 0.06; 1.2 0.02; 1.6 0.005];
L.uu.LS{3}  = [0 0.11; 0.2 0.13; 0.4 0.16; 0.6 0.17; 0.8 0.12; 1.0 0.07; 1.2 0.03; 1.6 0.01];
L.uu.DNS{4} = [0 0.045; 0.2 0.06; 0.35 0.085; 0.5 0.09; 0.7 0.075; 1.0 0.04; 1.5 0.01; 2.0 0.002];
L.uu.LS{4}  = [0 0.06; 0.2 0.07; 0.4 0.08; 0.6 0.088; 0.8 0.07; 1.0 0.05; 1.4 0.02; 1.6 0.01];

% ---- Figs. 3.25/3.26: <v'v'>/U0^2 profiles (upper half) ----
L.vv.DNS{1} = [0 0.12; 0.3 0.11; 0.5 0.07; 0.7 0.02; 1.0 0.005; 1.5 0];
L.vv.LS{1}  = [0 0.14; 0.3 0.15; 0.5 0.10; 0.7 0.04; 1.0 0.01; 1.5 0];
L.vv.DNS{2} = [0 0.45; 0.2 0.40; 0.4 0.26; 0.6 0.12; 0.8 0.05; 1.0 0.02; 1.5 0.002];
L.vv.LS{2}  = [0 0.35; 0.2 0.32; 0.4 0.22; 0.6 0.11; 0.8 0.05; 1.0 0.02; 1.5 0.005];
L.vv.DNS{3} = [0 0.52; 0.2 0.47; 0.4 0.33; 0.6 0.18; 0.8 0.08; 1.0 0.03; 1.5 0.005];
L.vv.LS{3}  = [0 0.32; 0.2 0.30; 0.4 0.23; 0.6 0.14; 0.8 0.07; 1.0 0.03; 1.5 0.005];
L.vv.DNS{4} = [0 0.36; 0.3 0.30; 0.6 0.19; 0.9 0.10; 1.2 0.05; 1.6 0.02; 2.0 0.005];
L.vv.LS{4}  = [0 0.23; 0.3 0.21; 0.6 0.15; 0.9 0.09; 1.2 0.05; 1.6 0.02];

% ---- Figs. 3.28/3.29: <u'v'>/U0^2 profiles (upper half; antisymmetric, negative for y>0) ----
L.uv.DNS{1} = [0 0; 0.2 0.02; 0.35 0.0; 0.45 -0.04; 0.55 -0.09; 0.62 -0.03; 0.75 0; 1.0 0];
L.uv.LS{1}  = [0 0; 0.2 0.01; 0.35 -0.01; 0.45 -0.04; 0.52 -0.08; 0.6 -0.04; 0.75 0; 1.0 0];
L.uv.DNS{2} = [0 0; 0.15 -0.04; 0.3 -0.09; 0.45 -0.13; 0.6 -0.08; 0.8 -0.02; 1.0 0];
L.uv.LS{2}  = [0 0; 0.15 -0.03; 0.3 -0.07; 0.45 -0.11; 0.6 -0.08; 0.8 -0.02; 1.0 0];
L.uv.DNS{3} = [0 0; 0.2 -0.06; 0.45 -0.12; 0.7 -0.07; 0.9 -0.03; 1.2 0];
L.uv.LS{3}  = [0 0; 0.2 -0.04; 0.45 -0.08; 0.7 -0.06; 0.9 -0.03; 1.2 -0.005];
L.uv.DNS{4} = [0 0; 0.2 -0.02; 0.6 -0.04; 0.9 -0.025; 1.2 -0.01; 1.6 0];
L.uv.LS{4}  = [0 0; 0.2 -0.02; 0.6 -0.045; 0.9 -0.03; 1.2 -0.01; 1.6 0];

% ---- qualitative statements used in the comparison ----
% Fig. 3.13: u_rms has two maxima per side (shear layer, strongest; second one slightly
%   downstream and closer to the axis = vortex formation); v_rms maximum on the axis
%   slightly downstream of bubble closure (x/D ~ 2.1 in the DNS).
L.vmax_x = 2.1;
% Sec. 3.4: centreline spectra -> v peaks at St (plus 3 St), u peaks at 2 St; -5/3 range.
% Sec. 4.5: U-/V-shape at x/D=1.06 depends on Lr; profiles collapse at equal x'/Lr,
%   x' = x - D/2; stations used x'/Lr = 0.431, 0.80, 1.17.  Averaging 50 vs 300 D/U0
%   changes Lr; intermittency in shear layer lasts ~10 cycles -> ~50 cycles needed.
% Sec. 3.3/6: plateau of centreline u at x/D~3 in experiments, absent in simulations;
%   hypothesised cause: vortex dislocations (Williamson), not resolved in pi*D span.
L.xLr = [0.431 0.80 1.17];
end
