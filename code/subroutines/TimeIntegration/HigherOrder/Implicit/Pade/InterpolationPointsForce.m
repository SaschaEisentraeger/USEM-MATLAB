function [xi] = InterpolationPointsForce(porder,Type)
%-------------------------------------------------------------------------%
% University of New South Wales (UNSW Sydney)
% Written by Sascha Eisentraeger (02/2021)
%-------------------------------------------------------------------------%
% Input:
% porder:           interpolation order of the force term
% Type:       		type of interpolation points; the boundaries of the
% 					interval are included in the first three options, while
% 					the fourth option features only internal points:
% 					(1) GLL-points;
% 					(2) GLC-points;
% 					(3) Babuska-Chen points;
%					(4) Chebyshev points;
%
% Output:
% xi:    			local coordinates of the interpolation points 
%               
%-------------------------------------------------------------------------%

switch Type
    case 1
        % GLL-points in the interval [-1,1]
        xi = lglnodes(porder);
    case 2
        % GLC-points in the interval [-1,1]
        xi = cgl_nodes(porder);
	case 3
        % Babuska-Chen points in the interval [-1,1]
        xi = bc_nodes(porder);
	case 4
		% Chebyshev points in the interval [-1,1]
		D = 2;
        xi = 1/2*D*(1+cos((2*(0:porder)'+1)./(2*porder+2)*pi)) - D/2;
end

xi = flip(xi);

end