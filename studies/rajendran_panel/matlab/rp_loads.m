function F = rp_loads(model, prob)
%RP_LOADS Equilibrated nodal forces on x = 0 and x = 2.
%   Consistent nodal forces of the exact tractions t = sigma n (n = +e_x on
%   x = 2, -e_x on x = 0), integrated exactly with the edge shape functions
%   of the spectral elements. For the 1 x 1, 2 x 2 and 4 x 4 meshes these
%   are the nodal forces and reactions listed in Rajendran (2010), Fig. 8
%   and Sec. 5.6 (checked in rp_verify).
p = model.p;
T = model.Topology_Node;
F = zeros(model.num_DoF, 1);
Qd = cyl_quad(p, 3*p.pxi + 4);
wq = Qd.wq;
for e = 1:size(T, 1)
    ids = T(e, :).';
    x = model.xCoor(ids);  y = model.yCoor(ids);
    for k = 1:4
        Ne = Qd.edge{k}.N;
        xe = Ne.'*x;  ye = Ne.'*y;
        for X = [0 2]
            if all(abs(xe - X) < 1e-12)
                if mod(k, 2) == 1
                    ds = abs(Qd.edge{k}.dNxi.'*y);
                else
                    ds = abs(Qd.edge{k}.dNeta.'*y);
                end
                sgn = 1 - 2*(X == 0);
                tx = sgn*prob.sxx(xe, ye);
                ty = sgn*prob.sxy(xe, ye);
                F(2*ids-1) = F(2*ids-1) + Ne*(tx.*wq.*ds);
                F(2*ids)   = F(2*ids)   + Ne*(ty.*wq.*ds);
            end
        end
    end
end
end
