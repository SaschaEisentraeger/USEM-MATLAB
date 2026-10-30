function [nodes,elements] = Read_mesh_ABAQUS(inp_file)
%READ_MESH_ABAQUS Nodes and 8-node elements of an ABAQUS input file.
%
%   [nodes, elements] = Read_mesh_ABAQUS(inp_file) reads the *Node block
%   and the first *Element block (type CPS8) of the file.
%
%   nodes:    [n x 2] x- and y-coordinates; the row number is the node label
%   elements: [m x 8] element topology

lines = splitlines(fileread(inp_file));
iNode = find(startsWith(lines,'*Node'),1);
iElem = find(startsWith(lines,'*Element'),1);
iEnd  = iElem + find(startsWith(lines(iElem+1:end),'*'),1);

nodes = sscanf(strjoin(lines(iNode+1:iElem-1),','),'%f,',[3 Inf])';
nodes = nodes(:,2:3);

elements = sscanf(strjoin(lines(iElem+1:iEnd-1),','),'%f,',[9 Inf])';
elements = elements(:,2:9);

end
