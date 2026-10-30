function M = ms_Me(P, x, y, nq, form)
%MS_ME Consistent element mass matrix (rho = 1) with nq x nq Gauss points:
%parametric test functions, trial functions of the formulation (form as in
%ms_Ke). Same integrand as ms_assemble_KM.
Qd = cyl_quad(P, nq);
[~, ~, ~, ~, detJ] = cyl_jacobian(Qd.dNxi, Qd.dNeta, x, y);
Np = Qd.N;
if strcmp(form, 'PP')
    Nt = Np;
else
    P.metricMethod = 1 + strcmp(form, 'PH');
    Nt = ShapeFunctions2D_Metric(Np.'*x, Np.'*y, x, y, P, [], [Qd.xi, Qd.eta]);
end
m = (Np.*(Qd.W(:).'.*detJ(:).'))*Nt.';
M = kron(m, eye(2));
end
