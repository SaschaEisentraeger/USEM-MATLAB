function results = run_tests(varargin)
%RUN_TESTS Run the USEM test suite.
%
%   results = run_tests()                   unit and integration tests
%   results = run_tests('unit')             only tests/unit
%   results = run_tests('integration')      only tests/integration
%   results = run_tests(..., 'KnownIssues', true)
%                                           also run tests tagged
%                                           'KnownIssue' (expected to fail
%                                           until the bugs listed in
%                                           docs/known-issues.md are fixed)
%   results = run_tests(..., 'Coverage', true)
%                                           write a Cobertura coverage
%                                           report for code/subroutines
%   results = run_tests(..., 'Parallel', true)
%                                           allow parfor loops in the code
%                                           to start a parallel pool
%
%   JUnit XML (and the coverage report) are written to tests/results/.
%   The function errors if any test fails, so it can be used in CI:
%       matlab -batch "cd tests; run_tests"

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.XMLPlugin
import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoberturaFormat
import matlab.unittest.selectors.HasTag

ip = inputParser;
ip.addOptional('scope', 'all', @(s) any(strcmpi(s, {'all', 'unit', 'integration'})));
ip.addParameter('KnownIssues', false, @islogical);
ip.addParameter('Coverage', false, @islogical);
ip.addParameter('Parallel', false, @islogical);
ip.parse(varargin{:});
opt = ip.Results;

testDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(testDir);
outDir  = fullfile(testDir, 'results');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% Several routines use parfor. Without an open pool MATLAB would start one
% automatically, which costs more time than the tests themselves. Disable
% the automatic pool creation for the duration of the run (restored below).
if ~opt.Parallel && ~isempty(ver('parallel'))
    ps = parallel.Settings;
    autoCreate = ps.Pool.AutoCreate;
    ps.Pool.AutoCreate = false;
    restorePool = onCleanup(@() set(ps.Pool, 'AutoCreate', autoCreate));
end

switch lower(opt.scope)
    case 'unit'
        suite = TestSuite.fromFolder(fullfile(testDir, 'unit'));
    case 'integration'
        suite = TestSuite.fromFolder(fullfile(testDir, 'integration'));
    otherwise
        suite = [TestSuite.fromFolder(fullfile(testDir, 'unit')), ...
                 TestSuite.fromFolder(fullfile(testDir, 'integration'))];
end
if ~opt.KnownIssues
    suite = suite.selectIf(~HasTag('KnownIssue'));
end

runner = TestRunner.withTextOutput;
runner.addPlugin(matlab.unittest.plugins.DiagnosticsRecordingPlugin);
runner.addPlugin(XMLPlugin.producingJUnitFormat(fullfile(outDir, 'junit.xml')));
if opt.Coverage
    runner.addPlugin(CodeCoveragePlugin.forFolder( ...
        fullfile(rootDir, 'code', 'subroutines'), 'IncludingSubfolders', true, ...
        'Producing', CoberturaFormat(fullfile(outDir, 'coverage.xml'))));
end

results = runner.run(suite);
failed = results([results.Failed]);
if ~isempty(failed)
    fprintf('\nFailed tests:\n');
    fprintf('  %s\n', failed.Name);
end
fprintf('\n%d passed, %d failed, %d incomplete (%.1f s)\n', ...
    nnz([results.Passed]), nnz([results.Failed]), ...
    nnz([results.Incomplete]), sum([results.Duration]));

% On GitHub Actions, report every failure as an annotation (visible on the
% run page without downloading the log)
if strcmp(getenv('GITHUB_ACTIONS'), 'true')
    githubAnnotations(failed);
    fprintf('::notice title=Test summary::%d passed, %d failed, %d incomplete\n', ...
        nnz([results.Passed]), nnz([results.Failed]), nnz([results.Incomplete]));
end

if any([results.Failed])
    error('run_tests:failed', '%d test(s) failed.', nnz([results.Failed]));
end

end

function githubAnnotations(failed)
% Report failures as GitHub workflow annotations. GitHub shows at most 10
% error annotations per step, so the names of all failed tests go into one
% annotation and the diagnostics of the first failure of each test class
% into the others.
if isempty(failed)
    return
end
esc = @(s) strrep(strrep(strrep(s, '%', '%25'), sprintf('\r'), ''), newline, '%0A');
names = {failed.Name};
fprintf('::error title=Failed tests::%s\n', esc(strjoin(names, newline)));
classes = regexprep(names, '/.*$', '');
[~, first] = unique(classes, 'stable');
for k = first(1:min(end, 9))'
    msg = '';
    rec = failed(k).Details.DiagnosticRecord;
    if ~isempty(rec)
        msg = rec(1).Report;
    end
    msg = msg(1:min(end, 4000));
    title = strrep(strrep(names{k}, ',', ' '), ':', ' ');
    fprintf('::error title=%s::%s\n', esc(title), esc(msg));
end
end
