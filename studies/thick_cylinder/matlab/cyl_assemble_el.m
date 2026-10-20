function Ke = cyl_assemble_el(model, C, form, e, nq)
%CYL_ASSEMBLE_EL Element stiffness matrix of element e (diagnostics).
p = model.p;
if nargin < 5 || isempty(nq), nq = 2*p.pxi; end
p.metricMethod = 1 + strcmp(form, 'PH');
ids = model.Topology_Node(e, :);
x = model.xCoor(ids);  y = model.yCoor(ids);
if any(strcmp(form, {'PM', 'PH'}))
    Ke = cyl_Ke_PM(p, C, x, y, nq);
else
    Ke = cyl_Ke_PP(p, C, x, y, nq);
end
end
