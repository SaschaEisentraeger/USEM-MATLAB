function R = ms_spectrum(pvec)
%MS_SPECTRUM Patch test, smooth field and generalized eigenvalues on the
%distorted 3 x 2 meshes (manuscript, Secs. 5.4-5.5). Wrapper around
%studies/metric_basis/study_hybrid_space.m; writes
%manuscript/data/spectrum-<mesh>.dat.
if nargin < 1, pvec = 1:10; end
[root, dataDir] = ms_setup();
addpath(fullfile(root, 'studies', 'metric_basis'));
R = study_hybrid_space('pvec', pvec, 'writeData', false);
meshes = fieldnames(R);
for k = 1:numel(meshes)
    ms_write_table(fullfile(dataDir, ['spectrum-' meshes{k} '.dat']), R.(meshes{k}));
end
end
