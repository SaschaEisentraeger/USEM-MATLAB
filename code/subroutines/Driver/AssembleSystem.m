function [K,M,R,dte,t_mat,ElemData] = AssembleSystem(model,settings,flag_comp)
%ASSEMBLESYSTEM Global system matrices of a model (wrapper of AssemblyFEM).
%
%   [K, M, R, dte, t_mat] = AssembleSystem(model, settings, flag_comp)
%   assembles the stiffness (and mass, damping) matrices of the model
%   returned by ReadModel, including the material data and Rayleigh
%   parameters added by RunSimulation (fields C, rho, zeta,
%   alpha_Rayleigh, beta_Rayleigh).
%
%   flag_comp: 0 K only; 1 K, M, R; 2 K, M, R and the critical time step
%              dte of the explicit methods (smallest over all elements)
%
%   [..., ElemData] = AssembleSystem(...) also returns the element
%   matrices (adaptive hybrid explicit method 109).

args = {model.Topology_Node,model.p,model.xCoor,model.yCoor,...
        model.num_Elem,model.num_DoF,model.C,model.rho,numel(model.C),...
        settings.mass_matrix,model.material_elset,model.alpha_Rayleigh,...
        model.beta_Rayleigh,model.zeta,flag_comp};

% AssemblyFEM computes the element data only if it is requested
if nargout > 5
    [K,M,R,dte,t_mat,ElemData] = AssemblyFEM(args{:});
else
    [K,M,R,dte,t_mat] = AssemblyFEM(args{:});
end

end
