classdef CodeHygieneTest < matlab.unittest.TestCase
    %CODEHYGIENETEST Conventions that keep the code base usable.
    %
    %   - Every file in code/ has help text ('help <name>' works).
    %   - No two files in code/ have the same name (the code folders are
    %     added to the path with genpath, so one would shadow the other).
    %   - Every entry of TimeIntegrationMethods points to a function of
    %     the code base.

    properties
        files
        codeDir
    end

    methods (TestClassSetup)
        function listFiles(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            tc.codeDir = fullfile(root, 'code');
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                tc.codeDir, 'IncludingSubfolders', true));
            tc.files = dir(fullfile(tc.codeDir, '**', '*.m'));
        end
    end

    methods (Test)
        function everyFileHasHelpText(tc)
            missing = {};
            for k = 1:numel(tc.files)
                f = fullfile(tc.files(k).folder, tc.files(k).name);
                if isempty(strtrim(help(f)))
                    missing{end+1} = f; %#ok<AGROW>
                end
            end
            tc.verifyEmpty(missing, 'Files without help text');
        end

        function fileNamesAreUnique(tc)
            names = {tc.files.name};
            [u, ~, j] = unique(lower(names));
            dup = u(accumarray(j(:), 1) > 1);
            tc.verifyEmpty(dup, 'Files that shadow each other');
        end

        function everySubroutineHasItsOwnFile(tc)
            % No local or nested functions: one function per file, so that
            % every subroutine can be found, called and tested on its own.
            % Third-party files are kept as published.
            thirdParty = {'quadmesh.m', 'expm3.m', 'maxdistcolor.m', ...
                          'maxdistcolor_view.m', 'srgb_to_Lab.m'};
            offenders = {};
            for k = 1:numel(tc.files)
                if any(strcmp(tc.files(k).name, thirdParty))
                    continue
                end
                txt = fileread(fullfile(tc.files(k).folder, tc.files(k).name));
                n = numel(regexp(txt, '^\s*function\s', 'lineanchors'));
                if n > 1
                    offenders{end+1} = tc.files(k).name; %#ok<AGROW>
                end
            end
            tc.verifyEmpty(offenders, 'Files with local or nested functions');
        end

        function timeIntegrationMethodsExist(tc)
            methods = TimeIntegrationMethods();
            for k = 1:numel(methods)
                f = which(func2str(methods(k).handle));
                tc.verifyTrue(startsWith(f, tc.codeDir), ...
                    sprintf('method %d: %s', methods(k).code, func2str(methods(k).handle)));
            end
        end
    end
end
