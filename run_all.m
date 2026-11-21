% RUN_ALL  Complete cylinder-wake PIV pipeline, in dependency order.
% 按依赖顺序运行完整的圆柱尾迹 PIV 处理流程。
%
% Before running / 运行前:
%   1. Put the course data in data/FM-Lab-WS2627-Group1/ (see data/README.md),
%      or call setenv('PIV_DATA', '<folder>') first.
%      把课程数据放到 data/FM-Lab-WS2627-Group1/，或先用 setenv 指定 PIV_DATA。
%   2. Optional: setenv('PIV_RESULTS', '<folder>') to write results elsewhere.
%      可选：用 PIV_RESULTS 指定结果目录。
%
% Runtime on a desktop PC: about 3 min (tested on R2025b).
% 普通台式机约 3 分钟（R2025b 上测试）。

repo = fileparts(mfilename('fullpath'));
addpath(repo, fullfile(repo, 'src'));

steps = { ...
    's1_calibration'      % magnification from the mean particle image (silhouette) / 由粒子平均图标定放大率
    's1b_ring'            % circle fit to the illuminated rim on the calibration image (D in px used by s2) / 标定图亮环拟合
    's2_process'          % displacement -> velocity, outlier filter, statistics / 位移转速度、剔除坏矢量、统计量
    's4_analysis'         % bulk flow, POD, phase averaging, spectra / 来流、POD、相位平均、频谱
    's5_figures'          % figures F01-F16 / 基本图
    's6_flowrate_timing'  % discharge Q, frame-rate check, F17-F18 / 流量与采样率检验
    's7_literature'       % comparison with Tremblay (2001) + new analyses / 文献对比与新分析
    's9_review_check'     % independent checks (needed by s8) / 复核计算（s8 需要）
    's8_lit_figures'      % figures L01-L03, N01-N06 / 对比图
    's10_fs_scenario'     % what if f_s = 2 Hz, N07 / 2 Hz 情景
    };

for k = 1:numel(steps)
    run_step(steps{k});
end
fprintf('\nAll steps finished. Results: %s\n', piv_path('results'));
