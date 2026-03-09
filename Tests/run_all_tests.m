% run_all_tests.m
% Run the full DeePhys unit test suite and print a pass/fail summary.
%
% Usage (from any directory):
%   run('/path/to/DeePhys/Tests/run_all_tests.m')
%
% Or after cd into the DeePhys root:
%   run('Tests/run_all_tests.m')
%
% Exit status: errors if any tests fail (non-zero exit code in batch mode).

root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'Functions')));
addpath(genpath(fullfile(root, 'Classes')));
addpath(genpath(fullfile(root, 'Toolboxes')));

test_dir = fileparts(mfilename('fullpath'));
suite    = matlab.unittest.TestSuite.fromFolder(test_dir, 'IncludingSubfolders', false);
runner   = matlab.unittest.TestRunner.withTextOutput('Verbosity', 2);
results  = runner.run(suite);

n_pass = sum([results.Passed]);
n_fail = sum([results.Failed]);
n_skip = sum([results.Incomplete]);

fprintf('\n%s\n', repmat('=', 1, 50));
fprintf('  DeePhys test suite results\n');
fprintf('%s\n', repmat('=', 1, 50));
fprintf('  Passed  : %i\n', n_pass);
fprintf('  Failed  : %i\n', n_fail);
fprintf('  Skipped : %i\n', n_skip);
fprintf('%s\n\n', repmat('=', 1, 50));

if n_fail > 0
    failing = results([results.Failed]);
    fprintf('Failing tests:\n');
    for k = 1:numel(failing)
        fprintf('  - %s\n', failing(k).Name);
    end
    fprintf('\n');
    error('DeePhys:testFailure', '%i test(s) failed.', n_fail);
end
