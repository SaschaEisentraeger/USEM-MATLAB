function ms_quadrature_calibration(family, pvec)
%MS_QUADRATURE_CALIBRATION Calibrated Gauss rules of the manuscript studies.
%
%   ms_quadrature_calibration(family, pvec) determines, for every
%   formulation (PP, PM, PH) and degree p in pvec, the number of Gauss
%   points per direction for which the element matrices of the mesh family
%   are integrated to round-off (ms_calibrate, tolerance 1e-12), and writes
%   manuscript/data/quadrature/<family>.dat with the columns
%       p, nG_<f>, conv_<f>, dev2p_<f>, rej_<f>, tol_<f>   (f = PP, PM, PH)
%   where conv = 0 marks rules that do not converge below the reference
%   rule max(64, 8p) (singular integrands on folded elements; nG is then the
%   reference rule), dev2p is the relative deviation of the former default
%   rule (2p points), rej = 1 marks trial spaces rejected by the code base
%   (singular interpolation matrix) and tol is the tolerance used
%   (ms_calibrate: 1e-12 or ten times the round-off level).
%
%   Families:
%     'conv010', 'conv020', 'conv030'  unit square, persistent curved-edge
%         distortion delta = 0.1, 0.2, 0.3 (ms_convergence); meshes with
%         N = 2, 4 and 8, maximum over the meshes (for N >= 4, the element
%         shapes repeat exactly); stiffness only
%     'curved15', 'curved20', 'combined30'  distorted 3 x 2 meshes
%         (ms_spectrum, ms_growth, ms_transient); stiffness and mass
%     'rect3x2'  undistorted 3 x 2 mesh (reference eigenvalues, PP only)
%     'cyl_aAAA_NN'  thick-cylinder mesh (alpha = AAA deg, N = NN); all
%         elements for N <= 2, otherwise the elements at the inner and outer
%         arc and the ring closest to r = 2; stiffness only

[root, dataDir] = ms_setup();
addpath(fullfile(root, 'studies', 'thick_cylinder', 'matlab'));
outDir = fullfile(dataDir, 'quadrature');
if ~exist(outDir, 'dir'), mkdir(outDir); end
forms = {'PP', 'PM', 'PH'};
withMass = false;
C = HookeMatrix2d(1, 0.3, 1);

if startsWith(family, 'conv')
    delta = str2double(family(5:end))/100;
    sets = {};
    for N = [2 4 8]
        [nodes, topo] = usemTestMesh(N, N, 1, 1, 2);
        nodes = usemDistortMesh(nodes, topo, 'curved', delta/N);
        sets{end+1} = @(p) meshElements(usemBuildModel(nodes, topo, p), []); %#ok<AGROW>
    end
elseif any(strcmp(family, {'curved15', 'curved20', 'combined30', 'rect3x2'}))
    withMass = true;
    [nodes, topo] = usemTestMesh(3, 2, 3, 2, 2);
    switch family
        case 'curved15',   nodes = usemDistortMesh(nodes, topo, 'curved', 0.15);
        case 'curved20',   nodes = usemDistortMesh(nodes, topo, 'curved', 0.20);
        case 'combined30', nodes = usemDistortMesh(nodes, topo, 'combined', 0.30);
    end
    if strcmp(family, 'rect3x2'), forms = {'PP'}; end
    sets = {@(p) meshElements(usemBuildModel(nodes, topo, p), [])};
elseif startsWith(family, 'cyl_')
    tok = sscanf(family, 'cyl_a%d_N%d');
    inp = fullfile(root, 'studies', 'thick_cylinder', 'mesh', 'inp', ...
                   sprintf('ThickCylinder_N%d_a%03d.inp', tok(2), tok(1)));
    E = 1000;  nu = 0.25;
    C = PlaneStressConstitutiveMatrices(HookeMatrix3d(E, nu));
    sets = {@(p) cylElements(cyl_model(inp, p), tok(2))};
else
    error('ms_quadrature_calibration:family', 'Unknown family "%s".', family);
end

T = struct('p', pvec(:));
for f = 1:numel(forms)
    for q = {'nG_', 'conv_', 'dev2p_', 'rej_', 'tol_'}
        T.([q{1} forms{f}]) = nan(numel(pvec), 1);
    end
end
for i = 1:numel(pvec)
    p = pvec(i);
    els = cellfun(@(s) s(p), sets, 'UniformOutput', false);
    for f = 1:numel(forms)
        t0 = tic;
        nG = 0;  conv = true;  d2p = 0;  rej = false;  tl = 0;
        for s = 1:numel(els)
            R = ms_calibrate(els{s}.xy, els{s}.P, C, forms{f}, withMass);
            if R.rejected, rej = true; continue, end
            nG = max(nG, R.nG);  conv = conv && R.conv;  d2p = max(d2p, R.dev2p);
            tl = max(tl, R.tol);
        end
        if rej && nG == 0, nG = NaN; conv = false; d2p = NaN; tl = NaN; end
        T.(['nG_' forms{f}])(i) = nG;
        T.(['conv_' forms{f}])(i) = conv;
        T.(['dev2p_' forms{f}])(i) = d2p;
        T.(['rej_' forms{f}])(i) = rej;
        T.(['tol_' forms{f}])(i) = tl;
        fprintf('%s p=%2d %s: nG = %3d  conv = %d  dev(2p) = %8.1e  tol = %8.1e  rej = %d  (%.0f s)\n', ...
                family, p, forms{f}, nG, conv, d2p, tl, rej, toc(t0));
    end
    ms_write_table(fullfile(outDir, [family '.dat']), T);
end
end

function S = meshElements(model, ids)
if isempty(ids), ids = 1:model.num_Elem; end
S.P = model.p;
S.xy = cell(numel(ids), 2);
for k = 1:numel(ids)
    nd = model.Topology_Node(ids(k), :);
    S.xy{k,1} = model.xCoor(nd);  S.xy{k,2} = model.yCoor(nd);
end
end

function S = cylElements(model, N)
if N <= 2
    ids = 1:model.num_Elem;
else
    % elements at the inner and outer arc and the ring closest to r = 2
    arc = unique(model.arcEdges(:, 1));
    rc = zeros(model.num_Elem, 1);
    for e = 1:model.num_Elem
        nd = model.Topology_Node(e, :);
        rc(e) = mean(hypot(model.xCoor(nd), model.yCoor(nd)));
    end
    [~, k] = sort(abs(rc - 2));
    ring = k(1:2*N);
    ids = unique([arc(:); ring(:)])';
end
S = meshElements(model, ids);
end
