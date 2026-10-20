function [K, M] = ms_assemble_KM(model, C, form)
%MS_ASSEMBLE_KM Dense global stiffness and consistent mass matrix for the
%classical (PP), Q_p(s,t) (PM) and hybrid (PH) spectral element.
P = model.p;
T = Topology_DoF_2D(model.Topology_Node);
n = model.num_DoF;
K = zeros(n); M = zeros(n);
[g, w] = lgwt(2*P.pxi, -1, 1);
[Xi, Eta] = meshgrid(g, g); [Wx, We] = meshgrid(w, w);
IP = [Xi(:), Eta(:)]; W = Wx(:).*We(:);
Np = ShapeFunctions2D_SEM(IP(:,1), IP(:,2), P);
for e = 1:model.num_Elem
    nd = model.Topology_Node(e, :);
    x = model.xCoor(nd); y = model.yCoor(nd);
    if strcmp(form, 'PP')
        Ke = usemKeSymmetric(P, C, x, y);
        Nt = Np;
    else
        P.metricMethod = 1 + strcmp(form, 'PH');
        Ke = Ke_2D(P, C, x, y);
        Xg = Local2GlobalCoordinatesQuad(P, IP, [x, y]);
        Nt = ShapeFunctions2D_Metric(Xg(:,1), Xg(:,2), x, y, P, [], IP);
    end
    Me = zeros(2*numel(nd));
    for k = 1:numel(W)
        J = Jacobi_Matrix(P, IP(k,:), x, y);
        Me = Me + kron(Np(:,k)*Nt(:,k)', eye(2))*det(J)*W(k);
    end
    K(T(e,:), T(e,:)) = K(T(e,:), T(e,:)) + Ke;
    M(T(e,:), T(e,:)) = M(T(e,:), T(e,:)) + Me;
end
end
