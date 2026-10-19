%% Compute stresses and strains for a 2D spectral element
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 27.11.2014 (University of Magdeburg)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% p:                Polynomial degree structure containing the field pxi
%                   and peta
% p.pxi:            Polynomial degree in xi
% p.peta:           Polynomial degree in eta
% xCoor:         	x-coordinates of the nodes of an an element
% yCoor:          	y-coordinates of the nodes of an an element
% u:                Element displacement vector
% C:                Elasticity matrix
% ndof:             Number of degrees of freedom per finite element
% local_coor:       Local coordinates of the points where the stresses and
% 					strains are computed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% stress:           Matrix of stresses
% strain:           Matrix of strains
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stress,strain] = Results_Stress_Strain(p,xCoor,yCoor,u,C,ndof,...
                           local_coor)

% Local coordinates
xi  = local_coor(:,1);
eta = local_coor(:,2);

% Number of points
NumP = numel(xi);

% Results
stress  = zeros(3,NumP);
strain  = zeros(3,NumP);

% % Derivatives of the parametric shape functions (xi,eta)
% [~,dNxi,dNeta] = ShapeFunctions2D_SEM(xi,eta,p);

% Derivatives of the metric shape functions (x,y)
IPg = Local2GlobalCoordinatesQuad(p,[xi,eta],[xCoor,yCoor]);
[~,dMx,dMy] = ShapeFunctions2D_Metric(IPg(:,1),IPg(:,2),xCoor,yCoor,p);

for i = 1:NumP
%------------------------Parametric shape functions-----------------------%
    % % B-Matrix: Strain-displacement matrix
    % % Jacobi matrix: Geometry mapping
    % [J] = Jacobi_Matrix(p,[xi(i) eta(i)],xCoor,yCoor);
    % % Derivative of the parametric shape function w.r.t. x and y
    % dN = J\[dNxi(:,i)';dNeta(:,i)'];
    % dNx = dN(1,:);
    % dNy = dN(2,:);
	% % Bu: parametric shape functions
    % Bu = zeros(3,ndof);
    % Bu(1,1:2:ndof) = dNx;
	% Bu(2,2:2:ndof) = dNy;
	% Bu(3,1:2:ndof) = dNy;
	% Bu(3,2:2:ndof) = dNx;
%--------------------------Metric shape functions-------------------------%
    % B-Matrix: Strain-displacement matrix
    % Bu: metric shape functions
	Bu = zeros(3,ndof);
	Bu(1,1:2:ndof) = dMx(:,i);
	Bu(2,2:2:ndof) = dMy(:,i);
	Bu(3,1:2:ndof) = dMy(:,i);
	Bu(3,2:2:ndof) = dMx(:,i);
%-------------------------------------------------------------------------%

    % Stress and strain values at the evaluation point
	strain(:,i) = Bu*u(1:ndof);
	stress(:,i) = C*strain(:,i);
end

end