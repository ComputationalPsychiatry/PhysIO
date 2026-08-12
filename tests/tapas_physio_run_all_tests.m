function testResults = tapas_physio_run_all_tests()
% Run the complete PhysIO unit and integration test suites.
%
%   testResults = tapas_physio_run_all_tests()
%
% Run tapas_physio_init() before calling this function so that PhysIO and
% SPM are initialized on the MATLAB path. If necessary, this function
% downloads both the PhysIO example data and the test reference results to
% their canonical PhysIO folders. Both suites are run recursively. An error
% is raised after both suites finish if any test failed or was incomplete,
% making this function suitable as a single pre-pull-request entry point
% for developers and coding agents.
%
% OUT
%   testResults    Struct with fields unit and integration containing the
%                  corresponding matlab.unittest.TestResult arrays.
%
% See also tapas_physio_run_unit_tests,
%          tapas_physio_run_integration_tests,
%          tapas_physio_get_paths_for_tests

% Individual test files do not perform this check when run directly; call
% tapas_physio_check_duplicate_paths first when bypassing an official runner.
tapas_physio_check_duplicate_paths();

% Verify both external datasets before starting either suite. If the
% example data or test reference results are missing, the helper downloads
% their version-matched Zenodo archives.
doUseZenodoPaths = true;
doVerifyPath = true;
doDownloadData = true;
[pathExamples, pathTestReferenceResults] = tapas_physio_get_paths_for_tests( ...
    doUseZenodoPaths, doVerifyPath, doDownloadData);

fprintf('PhysIO example data: %s\n', pathExamples);
fprintf('PhysIO test reference results: %s\n\n', pathTestReferenceResults);

fprintf('===== PhysIO unit tests =====\n');
testResults.unit = tapas_physio_run_unit_tests();

fprintf('\n===== PhysIO integration tests =====\n');
testResults.integration = tapas_physio_run_integration_tests();

fprintf('\n===== PhysIO test summary =====\n');
print_summary('Unit', testResults.unit);
print_summary('Integration', testResults.integration);

allResults = [testResults.unit, testResults.integration];
if isempty(allResults)
    error('tapas:physio:tests:NoTests', ...
        'No PhysIO unit or integration tests were discovered.');
end

if any([allResults.Failed]) || any([allResults.Incomplete])
    error('tapas:physio:tests:Failed', ...
        'PhysIO tests did not complete successfully. See the test summary above.');
end
end


function print_summary(label, results)
% Print compact totals for one test suite.

fprintf('%s: %d total, %d passed, %d failed, %d incomplete\n', ...
    label, numel(results), nnz([results.Passed]), ...
    nnz([results.Failed]), nnz([results.Incomplete]));
end
