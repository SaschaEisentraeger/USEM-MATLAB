function K = ms_Ke(P, C, x, y, nq, form)
%MS_KE Element stiffness matrix of the three formulations with nq x nq
%Gauss points.
%   form = 'PP' : conventional spectral element (cyl_Ke_PP, vectorized
%                 usemKeSymmetric)
%          'PM' : unsymmetric element, trial space Q_p in skew coordinates
%          'PH' : unsymmetric element, hybrid trial space
%   The unsymmetric matrices are computed with cyl_Ke_PM (vectorized
%   Ke_2D).
if strcmp(form, 'PP')
    K = cyl_Ke_PP(P, C, x, y, nq);
else
    P.metricMethod = 1 + strcmp(form, 'PH');
    K = cyl_Ke_PM(P, C, x, y, nq);
end
end
