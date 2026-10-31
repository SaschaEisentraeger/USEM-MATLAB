function [s, t, J0] = CorotationalSkewCoordinates(xQ, yQ, xN, yN, ~)
%COROTATIONALSKEWCOORDINATES Principal axes of the nodes instead of the skew
%coordinates. Takes the place of the routine of the code base for
%frame_cylinder_check only, which puts this folder first on the path.
xm = mean(xN);
ym = mean(yN);
[~, ~, W] = svd([xN(:) - xm, yN(:) - ym], 0);
Z = [xQ(:) - xm, yQ(:) - ym]*W;
s = Z(:, 1);
t = Z(:, 2);
J0 = W';
end
