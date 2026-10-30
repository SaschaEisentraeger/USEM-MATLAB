function [minDetJ, allDetJ, badMask] = MinJacobianDeterminant(X, T, type, minJac)
%MINJACOBIANDETERMINANT Smallest Jacobian determinant of all elements of a Q4/Q8 mesh.

if nargin<4, minJac = -inf; end
nE = size(T,1);
allDetJ = zeros(nE,1);
badMask = false(nE,1);
if type=="Q4"
    % 2x2 Gauss (±1/sqrt(3))
    g = 1/sqrt(3);
    gps = [-g -g; +g -g; +g +g; -g +g];
    for e=1:nE
        Xe = X(T(e,:),:);
        minj = inf;
        for q=1:4
            xi = gps(q,1); eta = gps(q,2);
            [dNdxi, dNdeta] = ShapeDerivativesQ4(xi,eta);
            J = [dNdxi'*Xe(:,1), dNdxi'*Xe(:,2); dNdeta'*Xe(:,1), dNdeta'*Xe(:,2)];
            dj = det(J);
            if dj < minj, minj = dj; end
        end
        allDetJ(e) = minj;
        badMask(e) = minj <= minJac;
    end
elseif type=="Q8"
    % 3x3 Gauss for safety
    g = [-sqrt(3/5), 0, +sqrt(3/5)];
    gps = combvec(g,g)'; % requires Neural Network TB in old MATLAB; roll our own:
    if size(gps,1)~=9
        gps = zeros(9,2); k=0; for a=1:3, for b=1:3, k=k+1; gps(k,:)=[g(a),g(b)]; end, end
    end
    for e=1:nE
        Xe = X(T(e,:),:);
        minj = inf;
        for q=1:size(gps,1)
            xi = gps(q,1); eta = gps(q,2);
            [dNdxi, dNdeta] = ShapeDerivativesQ8(xi,eta);
            J = [dNdxi'*Xe(:,1), dNdxi'*Xe(:,2); dNdeta'*Xe(:,1), dNdeta'*Xe(:,2)];
            dj = det(J);
            if dj < minj, minj = dj; end
        end
        allDetJ(e) = minj;
        badMask(e) = minj <= minJac;
    end
else
    error('Unsupported type');
end
minDetJ = min(allDetJ);
end
