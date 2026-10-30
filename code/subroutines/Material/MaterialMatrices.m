function [C,rho,zeta] = MaterialMatrices(material_properties,antype)
%MATERIALMATRICES Plane elasticity matrices of all materials of a model.
%
%   [C, rho, zeta] = MaterialMatrices(material_properties, antype) converts
%   the material data read by Read_material_ABAQUS into the input of
%   AssemblyFEM. Row i of material_properties describes material i:
%     {i,1} mass density, {i,3} elastic constants, {i,5} damping entry.
%   The elastic constants are either [E, nu] (isotropic) or the nine
%   entries of an orthotropic 6x6 stiffness matrix in the ABAQUS order
%   D1111 D1122 D2222 D1133 D2233 D3333 D1212 D1313 D2323
%   (Voigt order of this code: 11 22 33 12 23 13).
%
%   antype: 1 plane stress (sigma_33 = 0 condensed out), 2 plane strain.
%
%   C, rho, zeta: [nMat x 1] cell arrays with the 3x3 matrices, densities
%   and damping entries (empty if not given in the input file).

num_mat = size(material_properties,1);
C    = cell(num_mat,1);
rho  = cell(num_mat,1);
zeta = cell(num_mat,1);

for i = 1:num_mat
    constants = material_properties{i,3};
    if size(constants,2) == 2
        % Isotropic material
        C3d = HookeMatrix3d(constants(1),constants(2));
    else
        C3d = OrthotropicStiffness(constants);
    end

    switch antype
        case 1
            C{i} = PlaneStressConstitutiveMatrices(C3d);
        case 2
            C{i} = PlaneStrainConstitutiveMatrices(C3d);
        otherwise
            error('MaterialMatrices:antype', ...
                  'antype must be 1 (plane stress) or 2 (plane strain).');
    end
    rho{i}  = material_properties{i,1};
    zeta{i} = material_properties{i,5};
end

end
