function R = ms_growth(pvec, meshes, stepsPerPeriod)
%MS_GROWTH Growth rates of the free-vibration problem (manuscript, Sec. 5.7).
%
%   For M u'' + K u = 0, an eigenpair (lambda, phi) of (K, M) gives the
%   solutions exp(+-i sqrt(lambda) t) phi; any eigenvalue that is not real and
%   positive gives a growing solution with rate |Im sqrt(lambda)|. Written per
%   fundamental period T1 = 2 pi/omega_1 of the classical element:
%     gSemi = max |Im sqrt(lambda)| * T1          (growth exponent per period)
%   For the trapezoidal rule (Newmark 1/4, 1/2) with dt = T1/stepsPerPeriod the
%   amplification factors are z = (1 + i c)/(1 - i c), c = sqrt(lambda) dt/2
%   (both signs of the root):
%     gTrap = stepsPerPeriod * log(max |z|)        (growth exponent per period)
%   Mesh and boundary conditions as in ms_transient.
%   Writes manuscript/data/growth-<mesh>.dat.
if nargin < 1 || isempty(pvec), pvec = 2:8; end
if nargin < 2 || isempty(meshes), meshes = {'curved15', 'curved20', 'combined30'}; end
if nargin < 3 || isempty(stepsPerPeriod), stepsPerPeriod = 20; end
[~, dataDir] = ms_setup();
C = HookeMatrix2d(1, 0.3, 1);
forms = {'PP', 'PM', 'PH'};
for im = 1:numel(meshes)
    [nodes, topo] = usemTestMesh(3, 2, 3, 2, 2);
    switch meshes{im}
        case 'curved15',   nodes = usemDistortMesh(nodes, topo, 'curved', 0.15);
        case 'curved20',   nodes = usemDistortMesh(nodes, topo, 'curved', 0.20);
        case 'combined30', nodes = usemDistortMesh(nodes, topo, 'combined', 0.30);
    end
    T = struct('p', pvec(:));
    for f = 1:3
        T.(['gSemi_' forms{f}]) = nan(numel(pvec), 1);
        T.(['gTrap_' forms{f}]) = nan(numel(pvec), 1);
    end
    for i = 1:numel(pvec)
        model = usemBuildModel(nodes, topo, pvec(i));
        x = model.xCoor; y = model.yCoor;
        bnd = abs(x) < 1e-10 | abs(x - 3) < 1e-10 | abs(y) < 1e-10 | abs(y - 2) < 1e-10;
        act = setdiff((1:model.num_DoF)', sort([2*find(bnd) - 1; 2*find(bnd)]));
        for f = 1:3
            [K, M] = ms_assemble_KM(model, C, forms{f}, meshes{im});
            lam = eig(full(K(act, act)), full(M(act, act)));
            if f == 1
                T1 = 2*pi/sqrt(min(real(lam)));
                dt = T1/stepsPerPeriod;
            end
            w = sqrt(complex(lam));
            T.(['gSemi_' forms{f}])(i) = max(abs(imag(w)))*T1;
            c = [w; -w]*dt/2;
            z = abs((1 + 1i*c)./(1 - 1i*c));
            T.(['gTrap_' forms{f}])(i) = stepsPerPeriod*log(max(z));
        end
        fprintf('%-10s p=%d  gSemi: PP %8.2e PM %8.2e PH %8.2e | gTrap: PP %8.2e PM %8.2e PH %8.2e\n', ...
                meshes{im}, pvec(i), T.gSemi_PP(i), T.gSemi_PM(i), T.gSemi_PH(i), ...
                T.gTrap_PP(i), T.gTrap_PM(i), T.gTrap_PH(i));
    end
    ms_write_table(fullfile(dataDir, ['growth-' meshes{im} '.dat']), T);
end
end
