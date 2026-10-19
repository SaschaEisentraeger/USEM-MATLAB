function [K, M] = usemAssemble(model, C, rho, flag_comp, mass_matrix)
%USEMASSEMBLE Quiet wrapper around AssemblyFEM for a single material.
%
%   [K, M] = usemAssemble(model, C, rho, flag_comp, mass_matrix)
%   model       : struct from usemBuildModel
%   C           : 3x3 plane elasticity matrix
%   rho         : mass density (default 1)
%   flag_comp   : 0 = stiffness only (default), 1 = stiffness and mass
%   mass_matrix : 1 consistent (default), 2 nodal quadrature,
%                 3 row-sum, 4 HRZ

if nargin < 3 || isempty(rho),         rho = 1;         end
if nargin < 4 || isempty(flag_comp),   flag_comp = 0;   end
if nargin < 5 || isempty(mass_matrix), mass_matrix = 1; end

material_elset = {(1:model.num_Elem)'};

% AssemblyFEM prints progress information; evalc keeps the test log clean.
[~, K, M] = evalc(['AssemblyFEM(model.Topology_Node, model.p, ' ...
    'model.xCoor, model.yCoor, model.num_Elem, model.num_DoF, {C}, ' ...
    '{rho}, 1, mass_matrix, material_elset, 0, 0, {0}, flag_comp)']);

end
