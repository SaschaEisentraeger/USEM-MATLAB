function ms_nodes(p, shape)
%MS_NODES Node positions of a distorted element in physical, scaled skew
%and reference coordinates (manuscript, Fig. 1). Writes
%manuscript/data/nodes-<shape>-p<p>.dat (columns x y s t xi eta) and
%manuscript/data/outline-<shape>.dat (element boundary, columns x y s t).
if nargin < 1, p = 8; end
if nargin < 2, shape = 'general'; end
[~, dataDir] = ms_setup();
[x, y, P] = usemTestElement(shape, p);
[xi, eta] = ReferenceNodes2D(P);
[s, t] = CorotationalSkewCoordinates(x, y, x, y, P);
cs = (max(s) + min(s))/2; rs = (max(s) - min(s))/2;
ct = (max(t) + min(t))/2; rt = (max(t) - min(t))/2;
T = struct('x', x, 'y', y, 's', (s - cs)/rs, 't', (t - ct)/rt, 'xi', xi, 'eta', eta);
ms_write_table(fullfile(dataDir, sprintf('nodes-%s-p%d.dat', shape, p)), T);
% boundary of the element
r = linspace(-1, 1, 41)';
o = ones(size(r));
B = [r, -o; o, r; flipud(r), o; -o, flipud(r)];
Xb = Local2GlobalCoordinatesQuad(P, B, [x, y]);
[sb, tb] = CorotationalSkewCoordinates(Xb(:,1), Xb(:,2), x, y, P);
O = struct('x', Xb(:,1), 'y', Xb(:,2), 's', (sb - cs)/rs, 't', (tb - ct)/rt);
ms_write_table(fullfile(dataDir, sprintf('outline-%s.dat', shape)), O);
end
