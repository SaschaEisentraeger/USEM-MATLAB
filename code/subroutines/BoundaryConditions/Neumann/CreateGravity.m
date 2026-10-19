%% Volume load (2d) - Gravity load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 01/04/2019 (UNSW Sydney)
%
% Computation of the energetically equivalent force vector.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% p:                Polynomial degree structure containing the field pxi
%                   and peta
% p.pxi:            Polynomial degree in xi
% p.peta:           Polynomial degree in eta
% p.ndist:      	Nodal distribution
% x:                x-coordinates of all nodes
% y:                y-coordinates of all nodes
% num_DoF:          Total number of degrees of freedom
% TopoN:            Nodal topology (connectivity/coincidence) matrix
% grav:             Data of the gravity load
% mat:              Material data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% F_mech:           External load vector due to an applied gravity load
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F_mech] = CreateGravity(p,x,y,num_DoF,TopoN,grav,mat)

% Polynomial degrees in xi and eta
pxi = p.pxi;
peta = p.peta;

% Number of nodes per element
NumNodesElem = numel(TopoN(1,:));

% Number of elements in the model
nElem = size(TopoN,1);

% Element set
elset = grav.elset;
if isempty(elset)
    elset = (1:nElem);
end

% Number of loaded elements (the loops below run over the element set)
nElem = numel(elset);

% Number of different materials
NoMat = numel(mat.elset);

% Initialize the system load vector
F_mech = zeros(num_DoF,1);

% Elemental force vector
Fe = cell(nElem,1);

% Elemental degrees of freedom
DoFe = cell(nElem,1);

% Nodal coordinates (saved element-wise)
xCell = cell(nElem,1);
yCell = cell(nElem,1);

% Volume specific amplitude
amp_rho = grav.amplitude(2:3)*grav.amplitude(1);

% Mass density (saved element-wise)
rho_vec = zeros(nElem,1);

% Element data
for i = 1:nElem
    % Element number
    NoElem = elset(i);
    
    % Nodal coordinates
    xCell{i} = x(TopoN(NoElem,:));
    yCell{i} = y(TopoN(NoElem,:)); 
    
    % Degrees of freedom
    BC_aux = TopoN(NoElem,:);
    idx = zeros(1,2*NumNodesElem);
    idx(1:2:end) = BC_aux*2-1;
    idx(2:2:end) = BC_aux*2;
    DoFe{i} = idx;
    
    % Mass density for current element
    for m = 1:NoMat
    	if ~isempty(find(mat.elset{m} == NoElem,1))
            break;
     	end
    end
    rho_vec(i) = mat.prop{m,1};
end

% Determine the element force vectors due to the applied Neumann boundary
% conditions
parfor i = 1:nElem
    % Nodal coordinates
    xCoor = xCell{i};
    yCoor = yCell{i};
    
    % Volume specific amplitude
    amp_e = amp_rho*rho_vec(i);
    
    % Initialize the force vector
    Fb = zeros(2*NumNodesElem,1);
    
    % Loop over all integration points
    % Gauss-Legendre points (numerical integration)
    [xi,wxi]   = lgwt(pxi+1,-1,1);
    [eta,weta] = lgwt(peta+1,-1,1);

    % Save integration points in a vector
    [Xi,Eta] = meshgrid(xi,eta);
    IP = [reshape(Xi,[],1),reshape(Eta,[],1)];

    % Number of integration points
    NumIP = numel(IP(:,1));

    % Save integration weights in a vector
    [Wxi,Weta] = meshgrid(wxi,weta);
    IW = reshape(Wxi,[],1).*reshape(Weta,[],1);

    % Identity matrix
    I = eye(2);

    % Shape functions
    Naux = ShapeFunctions2D_SEM(IP(:,1),IP(:,2),p);
    
    for j = 1:NumIP
    	% Matrix of shape functions
        N = kron(Naux(:,j)',I);
        
        % Jacobi matrix: Geometry mapping
        [J] = Jacobi_Matrix(p,[IP(j,1) IP(j,2)],xCoor,yCoor);
        detJ = det(J);
        
        % Sum the contributions of each integration point to the element
        % force vector
        Fb = Fb + N'*amp_e*detJ*IW(j);
    end
    
    % Save the load
    Fe{i} = Fb;
end

% Assemble the system force vector
for i = 1:nElem
    F_mech(DoFe{i}) = F_mech(DoFe{i}) + Fe{i};
end

end