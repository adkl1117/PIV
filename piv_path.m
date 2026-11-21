function p = piv_path(which)
%PIV_PATH  Folder used by the analysis scripts, always with a trailing '/'.
%   piv_path('data')    -> folder that holds Calibration/, Data/ and Samples/
%   piv_path('results') -> folder where .mat results and figures/ are written
%
%   Defaults (relative to this repository):
%       data    = <repo>/data/FM-Lab-WS2627-Group1/
%       results = <repo>/results/
%   Override with environment variables, e.g. in MATLAB:
%       setenv('PIV_DATA',    'D:/my/PIV/FM-Lab-WS2627-Group1')
%       setenv('PIV_RESULTS', 'D:/my/PIV/results')
%
%   数据目录与结果目录。默认放在仓库内的 data/ 与 results/，
%   也可以用环境变量 PIV_DATA、PIV_RESULTS 指定其他位置。
%
%   Daten- und Ergebnisordner. Standardmäßig data/ und results/ im Repository;
%   andere Orte lassen sich über die Umgebungsvariablen PIV_DATA und PIV_RESULTS angeben.

repo = fileparts(mfilename('fullpath'));
switch lower(which)
    case 'data'
        p = getenv('PIV_DATA');
        if isempty(p), p = fullfile(repo, 'data', 'FM-Lab-WS2627-Group1'); end
        if ~isfolder(fullfile(p, 'Data'))
            error('piv_path:data', ['Data folder not found: %s\n' ...
                'Put the course data there (see data/README.md) or set PIV_DATA.'], p);
        end
    case 'results'
        p = getenv('PIV_RESULTS');
        if isempty(p), p = fullfile(repo, 'results'); end
        if ~isfolder(fullfile(p, 'figures')), mkdir(fullfile(p, 'figures')); end
    otherwise
        error('piv_path:arg', 'Use ''data'' or ''results''.');
end
p = strrep(p, '\', '/');
if p(end) ~= '/', p = [p '/']; end
end
