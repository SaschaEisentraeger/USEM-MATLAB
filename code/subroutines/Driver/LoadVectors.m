function [F,F_dt] = LoadVectors(model,flag_mRHS)
%LOADVECTORS Load vectors of a model (wrapper of NeumannBCs).
%
%   [F, F_dt] = LoadVectors(model) returns the constant part F and the
%   time-dependent part F_dt of the loads (concentrated forces, pressure,
%   traction, body forces, gravity) of the model returned by ReadModel.
%   F_dt is multiplied by the excitation signal s(t) in a transient
%   analysis.
%
%   [F, F_dt] = LoadVectors(model, flag_mRHS) with flag_mRHS = true returns
%   one column per load (multiple right-hand sides).

MAT.elset = model.material_elset;
MAT.prop  = model.material_properties;
args = {model.p,model.cload,model.dsload,model.body_force,model.gravity,...
        MAT,model.num_DoF,model.Topology_Node,model.xCoor,model.yCoor,...
        model.num_Nodes};
if nargin > 1
    args{end+1} = flag_mRHS;
end
[F,F_dt] = NeumannBCs(args{:});

end
