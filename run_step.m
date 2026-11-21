function run_step(name)
%RUN_STEP  Run one analysis script (src/<name>.m) in its own workspace.
%   Example: run_step('s6_flowrate_timing')
%   The scripts start with CLEAR, so each one runs inside a separate helper
%   function; the caller's variables (for example the loop in run_all.m) stay intact.
%   每个脚本开头都有 clear，放在单独的函数里运行，调用者的变量不会被清掉。
repo = fileparts(mfilename('fullpath'));
addpath(repo, fullfile(repo, 'src'));
t0 = tic;
fprintf('\n===== %s =====\n', name);
run_in_own_workspace(fullfile(repo, 'src', [name '.m']));
fprintf('----- %s done in %.0f s -----\n', name, toc(t0));
end

function run_in_own_workspace(script_file)
run(script_file);
end
