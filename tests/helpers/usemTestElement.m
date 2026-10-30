function [x, y, p] = usemTestElement(shape, pdeg)
%USEMTESTELEMENT Nodal coordinates of a single spectral element.
%
%   [x, y, p] = usemTestElement(shape, pdeg) returns the node
%   coordinates (in the order expected by ShapeFunctions2D_SEM) of one
%   element of degree pdeg, generated with pSpectralElement_NodesTopology.
%
%   shape = 'square'        : [0,2] x [0,1] rectangle
%           'parallelogram' : affine image of the reference square
%           'trapezoid'     : straight edges, angular distortion
%           'curved'        : 8-node input with two curved edges
%           'general'       : angular plus curved-edge distortion

switch shape
    case 'square'
        nodes = [0 0; 2 0; 2 1; 0 1];
    case 'parallelogram'
        nodes = [0 0; 2 0.3; 2.6 1.5; 0.6 1.2];
    case 'trapezoid'
        nodes = [0 0; 2 0; 1.4 1; 0.2 1.3];
    case 'curved'
        nodes = [0 0; 2 0; 2 1; 0 1; 1 -0.2; 2.25 0.5; 1 1; 0 0.5];
    case 'general'
        nodes = [0 0; 2 0.1; 1.5 1.2; 0.2 1.0; 1.1 -0.15; 1.9 0.6; ...
                 0.8 1.25; 0.05 0.5];
    otherwise
        error('usemTestElement:shape', 'Unknown element shape "%s".', shape);
end

topo = 1:size(nodes, 1);
model = usemBuildModel(nodes, topo, pdeg);
x = model.xCoor(model.Topology_Node(1,:));
y = model.yCoor(model.Topology_Node(1,:));
p = model.p;

end
