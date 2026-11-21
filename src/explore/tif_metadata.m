root = piv_path('data');
d = dir([root 'Samples/*.TIF']);
for k = [1 2 3 4]
  info = imfinfo([root 'Samples/' d(k).name]);
  fn = fieldnames(info);
  fprintf('--- %s\n', d(k).name);
  for j=1:numel(fn)
    v = info.(fn{j});
    if ischar(v) && ~isempty(v), fprintf('%s: %s\n', fn{j}, v(1:min(end,300))); end
  end
  if isfield(info,'UnknownTags'), for t=info.UnknownTags', fprintf('tag %d: %s\n', t.ID, mat2str(t.Value(1:min(end,40)))); end, end
end
