function cyl_verify(outfile)
%CYL_VERIFY Cross-checks of the study routines against the code base.
%
%   1. cyl_Ke_PM / cyl_Ke_PP vs Ke_2D / usemKeSymmetric (element level)
%   2. cyl_assemble vs AssemblyFEM (global stiffness)
%   3. Pressure load of the ABAQUS file (*Dsload, NeumannBCs/CreatePress)
%      vs the traction sigma_exact*n_h of cyl_loads (Lame problem)
%   4. Static solution with the code-base pipeline (AssemblyFEM, NeumannBCs,
%      DirichletBCs_homogeneous) vs cyl_assemble/cyl_solve with the same load
%   5. Geometry: signed area of the distorted meshes vs 2*pi
%   Writes a text report to outfile (default ../results/verification.txt).

here = fileparts(mfilename('fullpath'));
if nargin < 1, outfile = fullfile(here, '..', 'results', 'verification.txt'); end
cyl_setup_paths();
inpDir = fullfile(here, '..', 'mesh', 'inp');
inp = @(N, a) fullfile(inpDir, sprintf('ThickCylinder_N%d_a%03d.inp', N, a));
E = 1000;  nu = 0.25;
C = PlaneStressConstitutiveMatrices(HookeMatrix3d(E, nu));
lame = cyl_problem('lame', E, nu);
L = {};

L{end+1} = ['1. Element matrices: max_e ||K_study - K_code||_F/||K_code||_F, nq = 2p ' ...
            '(PM, PH: cyl_Ke_PM vs Ke_2D with p.metricMethod = 1, 2; PP: cyl_Ke_PP vs ' ...
            'usemKeSymmetric; "rejected": the code base stops on this mesh)'];
for a = [0 60 90]
    for pdeg = [2 5 8]
        m = cyl_model(inp(1, a), pdeg);
        r = [0 0 0];
        for e = 1:m.num_Elem
            ids = m.Topology_Node(e, :).';
            x = m.xCoor(ids);  y = m.yCoor(ids);
            for k = 1:2
                if isnan(r(k)), continue, end      % already rejected
                q = m.p;  q.metricMethod = k;
                try
                    K1 = Ke_2D(q, C, x, y, 2*pdeg, 2*pdeg);
                    r(k) = max(r(k), norm(cyl_Ke_PM(q, C, x, y, 2*pdeg) - K1, 'fro')/norm(K1, 'fro'));
                catch err
                    if ~startsWith(err.identifier, 'ShapeFunctions2D_Metric:'), rethrow(err); end
                    r(k) = NaN;
                end
            end
            K3 = usemKeSymmetric(m.p, C, x, y);
            r(3) = max(r(3), norm(cyl_Ke_PP(m.p, C, x, y, 2*pdeg) - K3, 'fro')/norm(K3, 'fro'));
        end
        c = arrayfun(@(v) sprintf('%.2e', v), r, 'UniformOutput', false);
        c(isnan(r)) = {'rejected'};
        L{end+1} = sprintf('   alpha = %2d, p = %d:  PM %s,  PH %s,  PP %s', a, pdeg, c{:}); %#ok<AGROW>
    end
end

L{end+1} = '2. Global stiffness: ||K_study - K_AssemblyFEM||_F/||K_AssemblyFEM||_F';
for a = [0 60]
    for pdeg = [2 4]
        m = cyl_model(inp(2, a), pdeg);
        K0 = codeStiffness(m, C);
        K = cyl_assemble(m, C, 'PM');
        L{end+1} = sprintf('   N = 2, alpha = %2d, p = %d:  %.2e', a, pdeg, ...
                           norm(full(K - K0), 'fro')/norm(full(K0), 'fro'));
    end
end

L{end+1} = ['3. Load vector (Lame): ||F_cyl_loads - F_NeumannBCs||/||F_NeumannBCs|| ' ...
            '(pressure of the ABAQUS file vs sigma_exact*n_h; differs by the ' ...
            'geometry error of the degree-p arc)'];
L{end+1} = '4. Static solution, same load (F_NeumannBCs): ||u_study - u_code||/||u_code||';
res3 = {};  res4 = {};
for a = [0 60]
    for pdeg = 2:8
        m = cyl_model(inp(1, a), pdeg);
        F0 = codeLoad(m);
        F = cyl_loads(m, lame);
        res3{end+1} = sprintf('   N = 1, alpha = %2d, p = %d:  %.2e', a, pdeg, norm(F - F0)/norm(F0)); %#ok<AGROW>
        % code-base pipeline (as in call_Static.m)
        K0 = codeStiffness(m, C);
        act0 = setdiff((1:m.num_DoF).', codeConstrainedDofs(m.boundary_conditions));
        u0 = zeros(m.num_DoF, 1);
        u0(act0) = K0(act0, act0)\F0(act0);
        % study routines with the code-base load vector
        try
            K = cyl_assemble(m, C, 'PM');
        catch err
            if ~startsWith(err.identifier, 'ShapeFunctions2D_Metric:'), rethrow(err); end
            res4{end+1} = sprintf('   N = 1, alpha = %2d, p = %d:  rejected by the code base', a, pdeg); %#ok<AGROW>
            continue
        end
        [~, ~, act] = cyl_solve(m, K, F, lame);
        u = zeros(m.num_DoF, 1);
        u(act) = K(act, act)\F0(act);
        res4{end+1} = sprintf('   N = 1, alpha = %2d, p = %d:  %.2e (same constrained DOFs: %d)', ...
                              a, pdeg, norm(u - u0)/norm(u0), isequal(sort(act(:)), sort(act0(:)))); %#ok<AGROW>
    end
end
i3 = find(strncmp(L, '3.', 2));
L = [L(1:i3), res3, L(i3+1:end), res4];

L{end+1} = '5. Geometry: signed area - 2*pi (p = 8), number of elements with det J <= 0';
for N = [1 2 4 8]
    for a = [0 15 30 45 60 90]
        m = cyl_model(inp(N, a), 8);
        L{end+1} = sprintf('   N = %d, alpha = %2d:  %+.2e,  %3d of %3d elements, min det J = %+.3e', ...
                           N, a, m.area - 2*pi, m.negDetJ, m.num_Elem, m.minDetJ); %#ok<AGROW>
    end
end

fid = fopen(outfile, 'w');
fprintf(fid, '%s\n', L{:});
fclose(fid);
fprintf('%s\n', L{:});
end

function K0 = codeStiffness(m, C)
% AssemblyFEM (code base): stiffness only, default quadrature of Ke_2D
mat = {(1:m.num_Elem).'};
evalc(['K0 = AssemblyFEM(m.Topology_Node, m.p, m.xCoor, m.yCoor, m.num_Elem, ' ...
       'm.num_DoF, {C}, {1}, 1, 1, mat, 0, 0, {0}, 0);']);
end

function dofs = codeConstrainedDofs(bc)
% Constrained DOFs from the ABAQUS boundary conditions (code base)
dofs = unique(DirichletBCs_homogeneous(bc));
end

function F0 = codeLoad(m)
% NeumannBCs (code base) with the *Dsload pressure of the ABAQUS file
MAT.elset = m.material_elset;
MAT.prop = m.material_properties;
evalc(['F0 = NeumannBCs(m.p, m.cload, m.dsload, m.body_force, m.gravity, MAT, ' ...
       'm.num_DoF, m.Topology_Node, m.xCoor, m.yCoor, m.num_Nodes);']);
end
