function make_abaqus_inputs()
%MAKE_ABAQUS_INPUTS Write the example ABAQUS input files.
%
%   make_abaqus_inputs() reconstructs CAE-style input files for the Lee &
%   Bathe patch tests in code/InputFilesAbaqus from the reference
%   workspaces in tests/fixtures/reference (mesh, material, and the
%   original node/element sets of the loads and boundary conditions). Read
%   with the ABAQUS readers of the code, the files reproduce the stored
%   model data exactly (see tests/integration/AbaqusInputTest.m).

here = fileparts(mfilename('fullpath'));
d   = [fullfile(here, 'reference') filesep];
out = [fullfile(fileparts(fileparts(here)), 'code', 'InputFilesAbaqus') filesep];
if ~exist(out, 'dir')
    mkdir(out);
end
cases = {'PatchTestLeeQuadratic_undistorted', 'ndof12', 'quadratic patch test (pure bending), undistorted mesh'; ...
         'PatchTestLeeQuadratic_angular',     'ndof12', 'quadratic patch test (pure bending), angular distortion'; ...
         'PatchTestLeeQuadratic_midedge',     'ndof12', 'quadratic patch test (pure bending), mid-edge (curved edge) distortion'; ...
         'PatchTestLeeCubic_undistorted',     'ndof30', 'cubic patch test, undistorted mesh'; ...
         'PatchTestLeeCubic_angular',         'ndof56', 'cubic patch test, angular distortion'; ...
         'PatchTestLeeCubic_midedge',         'ndof30', 'cubic patch test, mid-edge (curved edge) distortion'; ...
         'PatchTestLeeCubic_Rajendran_adj',   'ndof104','cubic patch test, five distorted 8-node elements (Rajendran & Liew 2003)'};
for c = 1:size(cases,1)
    S = load([d cases{c,1} '_' cases{c,2} '.mat']);
    X = S.nodes; T = S.topology; nn = size(X,1);
    if size(T,2) == 4, et = 'CPS4'; else, et = 'CPS8'; end
    % map an adapted node set back to the original mesh (by coordinates)
    origSet = @(ns) mapSet(ns, S, X);
    fid = fopen([out cases{c,1} '.inp'], 'w');
    fprintf(fid, '*Heading\n** Job name: %s Model name: Model-1\n', cases{c,1});
    fprintf(fid, '** Lee & Bathe (1993) %s.\n', cases{c,3});
    fprintf(fid, '** Plate [0,10] x [-1,1], E = 70 GPa, nu = 0.3, plane stress.\n');
    fprintf(fid, '** Spatially varying edge loads (NU) are defined in SpatialLoadingFunction.m.\n');
    fprintf(fid, '** Reconstructed from the result workspace %s_%s.mat.\n', cases{c,1}, cases{c,2});
    fprintf(fid, '*Preprint, echo=NO, model=NO, history=NO, contact=NO\n**\n** PARTS\n**\n*Part, name=Plate\n*Node\n');
    for i = 1:nn
        fprintf(fid, '%7d, %12.6g, %12.6g\n', i, X(i,1), X(i,2));
    end
    fprintf(fid, '*Element, type=%s\n', et);
    for e = 1:size(T,1)
        fprintf(fid, '%d', e); fprintf(fid, ', %d', T(e,:)); fprintf(fid, '\n');
    end
    fprintf(fid, '*Nset, nset=Set-all, generate\n 1, %d, 1\n', nn);
    fprintf(fid, '*Elset, elset=Set-all, generate\n 1, %d, 1\n', size(T,1));
    fprintf(fid, '** Section: Section-1\n*Solid Section, elset=Set-all, material=Aluminium\n1.,\n*End Part\n**\n**\n');
    fprintf(fid, '** ASSEMBLY\n**\n*Assembly, name=Assembly\n**\n*Instance, name=Plate-1, part=Plate\n*End Instance\n**\n');
    % boundary condition sets (stored in reverse order)
    nbc = numel(S.boundary_conditions);
    for k = nbc:-1:1
        bc = S.boundary_conditions(k);
        writeSet(fid, sprintf('Set-%d', nbc-k+1), origSet(bc.nset), bc.elset);
    end
    for k = 1:numel(S.dsload)
        q = S.dsload(k);
        writeSet(fid, regexprep(q.set_name, '\d+$', ''), origSet(q.nset), q.elset);
    end
    fprintf(fid, '*End Assembly\n**\n** MATERIALS\n**\n*Material, name=Aluminium\n*Density\n%g.,\n*Elastic\n %g, %g\n', ...
        S.material_properties{1,1}, S.material_properties{1,3}(1), S.material_properties{1,3}(2));
    fprintf(fid, '** ----------------------------------------------------------------\n**\n');
    fprintf(fid, '** STEP: Step-1\n**\n*Step, name=Step-1, nlgeom=NO\n*Static\n1., 1., 1e-05, 1.\n**\n** BOUNDARY CONDITIONS\n**\n');
    for k = nbc:-1:1
        bc = S.boundary_conditions(k);
        fprintf(fid, '** Name: BC-%d Type: Displacement/Rotation\n*Boundary\n', nbc-k+1);
        switch bc.type
            case 'fixed', fprintf(fid, 'Set-%d, ENCASTRE\n', nbc-k+1);
            case 'Ux 0',  fprintf(fid, 'Set-%d, 1, 1\n', nbc-k+1);
            case 'Uy 0',  fprintf(fid, 'Set-%d, 2, 2\n', nbc-k+1);
            otherwise, error('unsupported BC %s', bc.type);
        end
    end
    fprintf(fid, '**\n** LOADS\n**\n');
    for k = 1:numel(S.dsload)
        q = S.dsload(k);
        nm = regexprep(q.set_name, '\d+$', '');
        switch q.type
            case 'Pressure fx'
                fprintf(fid, '** Name: Load-%d   Type: Pressure\n*Dsload\n%s, PNU, %g.\n', k, nm, q.amplitude);
            case 'Traction fx'
                v = q.DirectionVector;
                fprintf(fid, '** Name: Load-%d   Type: Surface traction\n*Dsload\n%s, TRVECNU, %g., %g., %g., %g.\n', ...
                    k, nm, q.amplitude, v(1), v(2), v(3));
            otherwise, error('unsupported load %s', q.type);
        end
    end
    fprintf(fid, '**\n** OUTPUT REQUESTS\n**\n*Restart, write, frequency=0\n*Output, field, variable=PRESELECT\n*Output, history, variable=PRESELECT\n*End Step\n');
    fclose(fid);
    fprintf('%s written\n', cases{c,1});
end

end

function ns = mapSet(ns, S, X)
% coordinates of the (adapted) set nodes
P = [S.xCoor(ns) S.yCoor(ns)];
if numel(ns) == 1
    ns = find(abs(X(:,1)-P(1)) < 1e-9 & abs(X(:,2)-P(2)) < 1e-9);
    return
end
% all set nodes lie on one straight edge x = const
assert(all(abs(P(:,1) - P(1,1)) < 1e-9), 'set is not on a vertical edge');
ns = find(abs(X(:,1) - P(1,1)) < 1e-9);
end

function writeSet(fid, name, nset, elset)
fprintf(fid, '*Nset, nset=%s, instance=Plate-1\n', name);
fprintf(fid, ' %s\n', strjoin(compose('%d', double(nset(:))'), ', '));
if ~isempty(elset)
    fprintf(fid, '*Elset, elset=%s, instance=Plate-1\n', name);
    fprintf(fid, ' %d,', elset); fprintf(fid, '\n');
end
end
