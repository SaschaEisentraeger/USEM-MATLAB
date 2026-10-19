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
            tc.verifyEqual(Cps, HookeMatrix2d(tc.E, nu, 1), 'RelTol', 1e-12);
            ref = tc.E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
            tc.verifyEqual(Cps, ref, 'AbsTol', 1e-12*tc.E);
        end

        function planeStrainReductionMatchesClosedForm(tc, nu)
            Cpe = PlaneStrainConstitutiveMatrices(HookeMatrix3d(tc.E, nu));
            tc.verifyEqual(Cpe, HookeMatrix2d(tc.E, nu, 2), 'RelTol', 1e-12);
        end

        function hookeMatrix2dRejectsInvalidAnalysisType(tc)
            tc.verifyError(@() HookeMatrix2d(tc.E, 0.3, 3), ?MException);
        end

        function transverselyIsotropicReducesToIsotropic(tc)
            nu0 = 0.3; G = tc.E/(2*(1+nu0));
            C = HookeMatrix3d_transiso(tc.E, tc.E, G, nu0, nu0);
            tc.verifyEqual(C, HookeMatrix3d(tc.E, nu0), 'RelTol', 1e-10, ...
                'AbsTol', 1e-6*tc.E);
        end

        function rotationByZeroIsIdentity(tc)
            C = HookeMatrix3d_transiso(150e9, 10e9, 5e9, 0.3, 0.4);
            for ax = {'x', 'y', 'z'}
                tc.verifyEqual(RotationHooke(C, 0, ax{1}), C);
            end
        end

        function rotationPreservesIsotropy(tc)
            C = HookeMatrix3d(tc.E, 0.3);
            for ax = {'x', 'y', 'z'}
                tc.verifyEqual(RotationHooke(C, 0.6, ax{1}), C, ...
                    'AbsTol', 1e-10*tc.E, ['axis ' ax{1}]);
            end
        end

        function rotationForthAndBackIsIdentity(tc)
            C = HookeMatrix3d_transiso(150e9, 10e9, 5e9, 0.3, 0.4);
            for ax = {'x', 'y', 'z'}
                Cr = RotationHooke(RotationHooke(C, 0.4, ax{1}), -0.4, ax{1});
                tc.verifyEqual(Cr, C, 'AbsTol', 1e-10*150e9, ['axis ' ax{1}]);
            end
        end

        function transverselyIsotropicIsInvariantAboutFibreAxis(tc)
            % Fibre direction = 1 (x), isotropy plane = 2-3
            C = HookeMatrix3d_transiso(150e9, 10e9, 5e9, 0.3, 0.4);
            tc.verifyEqual(RotationHooke(C, 0.9, 'x'), C, 'AbsTol', 1e-10*150e9);
            tc.verifyGreaterThan(norm(RotationHooke(C, 0.9, 'z') - C), 1e9);
        end
    end
end
