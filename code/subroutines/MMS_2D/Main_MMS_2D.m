% Add to search path
disp('Add path')
path1 = genpath('.\code');
path2 = genpath('.\differentialOperators');
path3 = genpath('.\examples');
addpath(path1,path2,path3);

%% Call example
% Hauke:
% examples_cdrEquation
% examples_PoissonEquation
% examples_ReynoldsEquation
% examples_twoPhaseFlow
% examples_twoPhaseFlow_transient

% Sascha:
examples_Elasticity
% examples_Elasticity_transient