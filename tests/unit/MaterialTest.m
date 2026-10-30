classdef MaterialTest < matlab.unittest.TestCase
    %MATERIALTEST Constitutive matrices (Voigt notation).

    properties (TestParameter)
        nu = {0, 0.25, 0.3, 0.49};
    end

    properties (Constant)
        E = 210e9;
    end

    methods (TestClassSetup)
        function addCodeToPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'code', 'subroutines'), 'IncludingSubfolders', true));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'tests', 'helpers')));
        end
    end

    methods (Test)
        function isotropicHookeMatrixIsSymmetricPositiveDefinite(tc, nu)
            C = HookeMatrix3d(tc.E, nu);
            tc.verifySize(C, [6 6]);
            tc.verifyEqual(C, C', 'AbsTol', 0);
            tc.verifyGreaterThan(min(eig(C)), 0);
        end

        function isotropicHookeMatrixHasLameConstants(tc, nu)
            C = HookeMatrix3d(tc.E, nu);
            lambda = tc.E*nu/((1+nu)*(1-2*nu));
            mu = tc.E/(2*(1+nu));
            tc.verifyEqual(C(1,1), lambda + 2*mu, 'RelTol', 1e-14);
            tc.verifyEqual(C(1,2), lambda, 'AbsTol', 1e-14*tc.E);
            tc.verifyEqual(C(4,4), mu, 'RelTol', 1e-14);
        end

        function uniaxialStressGivesYoungsModulus(tc, nu)
            % Compliance S = inv(C): S(1,1) = 1/E, S(1,2) = -nu/E
            S = inv(HookeMatrix3d(tc.E, nu));
            tc.verifyEqual(S(1,1), 1/tc.E, 'RelTol', 1e-12);
            tc.verifyEqual(S(1,2), -nu/tc.E, 'AbsTol', 1e-12/tc.E);
        end

        function planeStressReductionMatchesClosedForm(tc, nu)
            Cps = PlaneStressConstitutiveMatrices(HookeMatrix3d(tc.E, nu));
            tc.verifyEqual(Cps, HookeMatrix2d(tc.E, nu), 'RelTol', 1e-12);
            ref = tc.E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
            tc.verifyEqual(Cps, ref, 'AbsTol', 1e-12*tc.E);
        end

    end
end
