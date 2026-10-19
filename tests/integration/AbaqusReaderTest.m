classdef AbaqusReaderTest < matlab.unittest.TestCase
    %ABAQUSREADERTEST Robustness of the ABAQUS input file readers.
    %
    %   Variants of code/InputFilesAbaqus/PatchTestLeeQuadratic_angular.inp
    %   are written to a temporary folder and read with the steps of
    %   Start_Sim.m.

    properties
        base        % content of the reference input file
    end

    methods (TestClassSetup)
        function addCodeToPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'code', 'subroutines'), 'IncludingSubfolders', true));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'tests', 'helpers')));
            tc.base = fileread(fullfile(root, 'code', 'InputFilesAbaqus', ...
                'PatchTestLeeQuadratic_angular.inp'));
        end
    end

    methods (TestMethodSetup)
        function useTemporaryFolder(tc)
            tc.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
        end
    end

    methods (Test)
        function referenceFileIsRead(tc)
            S = read(tc, tc.base);
            verifyLeftEdgeConstrained(tc, S);
        end

        function setNamesAreMatchedExactly(tc)
            % Set-10 (defined first) must not be taken for Set-1
            txt = edit(tc, tc.base, '*Nset, nset=Set-1, instance=Plate-1', ...
                sprintf(['*Nset, nset=Set-10, instance=Plate-1\n 1, 2\n' ...
                         '*Elset, elset=Set-10, instance=Plate-1\n 1,\n' ...
                         '*Nset, nset=Set-1, instance=Plate-1']));
            S = read(tc, txt);
            verifyLeftEdgeConstrained(tc, S);
        end

        function missingBoundaryConditionSetRaisesError(tc)
            % used to be an endless loop
            txt = edit(tc, tc.base, sprintf('*Boundary\nSet-2, 2, 2'), ...
                sprintf('*Boundary\nSet-7, 2, 2'));
            write(tc, txt);
            tc.verifyError(@() Read_boundary_condition_ABAQUS('model_adjusted.inp'), ...
                'Read_boundary_condition_ABAQUS:setNotFound');
        end

        function boundaryConditionOnTwoDirectionsInSeparateLines(tc)
            % the U2 branch used to loop without reading the next line
            txt = edit(tc, tc.base, sprintf('*Boundary\nSet-2, 2, 2'), ...
                sprintf('*Boundary\nSet-2, 2, 2\nSet-2, 3, 3'));
            S = read(tc, txt);
            tc.verifyEqual(S.boundary_conditions(1).type, 'UyUz 0');
        end

        function percentSignsAreCopiedLiterally(tc)
            txt = edit(tc, tc.base, '*Heading', ...
                sprintf('*Heading\n** 100%% load, path C:\\temp\\new'));
            write(tc, txt);
            adj = fileread('model_adjusted.inp');
            tc.verifySubstring(adj, '** 100% load, path C:\temp\new');
        end

        function nonConsecutiveNodeLabelsRaiseError(tc)
            txt = edit(tc, tc.base, sprintf('\n      6,'), sprintf('\n      7,'));
            write(tc, txt);
            tc.verifyError(@() Read_mesh_ABAQUS('model_adjusted.inp'), ...
                'Read_mesh_ABAQUS:nodeLabels');
        end

        function gravityWithThreeComponents(tc)
            txt = edit(tc, tc.base, sprintf('Pressure, PNU, 1.\n'), ...
                sprintf(['Pressure, PNU, 1.\n** Name: Load-2   Type: Gravity\n' ...
                         '*Dload\n, GRAV, 9.81, 0., -1., 0.\n']));
            S = read(tc, txt);
            tc.verifyEqual(S.gravity.amplitude, [9.81; 0; -1]);
        end

        function concentratedLoadOnGeneratedSet(tc)
            txt = edit(tc, tc.base, '*End Assembly', ...
                sprintf('*Nset, nset=Set-3, instance=Plate-1, generate\n 3, 4, 1\n*End Assembly'));
            txt = edit(tc, txt, sprintf('Pressure, PNU, 1.\n'), ...
                sprintf('Pressure, PNU, 1.\n** Name: Load-2   Type: Concentrated force\n*Cload\nSet-3, 2, -5.\n'));
            S = read(tc, txt);
            tc.verifyEqual(double(S.cload.nset), [3; 4]);
            tc.verifyEqual(S.cload.type, 'f_y');
            tc.verifyEqual(S.cload.amplitude, -5);
        end
    end
end

function verifyLeftEdgeConstrained(tc, S)
% Set-1 (u_x = 0): element 2 and all nodes of the left edge x = 0
% (corner nodes 5, 6 and the mid-side node of the p = 2 mesh)
bc = S.boundary_conditions(2);
tc.verifyEqual(bc.type, 'Ux 0');
tc.verifyEqual(double(bc.elset), 2);
tc.verifyEqual(S.xCoor(bc.nset), zeros(3, 1));
tc.verifyEqual(sort(S.yCoor(bc.nset)), [-1; 0; 1]);
end

function txt = edit(tc, txt, old, new)
% Replace exactly one occurrence (fail loudly if the base file changes)
txt = strrep(txt, sprintf('\r\n'), newline);
tc.assertEqual(numel(strfind(txt, old)), 1, ['pattern not found: ' old]);
txt = strrep(txt, old, new);
end

function write(~, txt)
fid = fopen('model.inp', 'w');
fprintf(fid, '%s', txt);
fclose(fid);
Adjust_Set_names_input_file_ABAQUS('model.inp');
end

function S = read(~, txt)
fid = fopen('model.inp', 'w');
fprintf(fid, '%s', txt);
fclose(fid);
p.pxi = 2; p.peta = 2; p.ndist = 2;
[~, S] = evalc('usemReadAbaqus(''model.inp'', p)');
end
