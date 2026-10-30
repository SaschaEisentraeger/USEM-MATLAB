function p = AlignPolynomialDegrees(Topology_Node,p,verbose)
%ALIGNPOLYNOMIALDEGREES Use pxi = peta if the element axes are not aligned.
%
%   p = AlignPolynomialDegrees(Topology_Node, p) checks every edge shared by
%   two elements. Edges 1 and 3 of an element are parallel to its xi axis,
%   edges 2 and 4 to its eta axis. If a shared edge is an xi edge of one
%   element and an eta edge of the other, the degrees pxi and peta of the
%   two elements would not match along that edge; both are then set to
%   min(pxi, peta). Nothing changes if p.pxi == p.peta.
%
%   p = AlignPolynomialDegrees(..., verbose) prints a message when the
%   degrees are changed (default true).

if nargin < 3
    verbose = true;
end

if p.pxi == p.peta
    return
end

[Topology_Edge,NumEdges] = Topology_Edge_2D(Topology_Node);
[EdgeElemList,LocalEdgeNum] = ConnectivityEdgeElement(Topology_Edge);
for i = 1:NumEdges
    if numel(EdgeElemList{i}) == 2
        % Odd local edge numbers belong to xi, even ones to eta
        r1 = mod(LocalEdgeNum{i}(1),2);
        r2 = mod(LocalEdgeNum{i}(2),2);
        if r1 ~= r2
            psel = min([p.pxi,p.peta]);
            p.pxi  = psel;
            p.peta = psel;
            if verbose
                disp([newline 'Local coordinate systems of adjacent' ...
                      ' elements are not aligned: pxi=peta=' ...
                      num2str(psel) newline])
            end
            return
        end
    end
end

end
