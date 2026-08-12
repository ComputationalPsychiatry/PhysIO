function tests = tapas_physio_print_figs_to_file_test()
% Tests figure-report export formats.

tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testFolder = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(testFolder)));
testCase.TestData.pathFixture = ...
    matlab.unittest.fixtures.PathFixture( ...
    fullfile(repoRoot, 'code'), IncludingSubfolders=true);
testCase.applyFixture(testCase.TestData.pathFixture);
end

function test_exports_multipage_pdf(testCase)
[verbose, testFolder] = create_test_figures(testCase, 'report.pdf');

verbose = tapas_physio_print_figs_to_file(verbose, testFolder);

verifyEqual(testCase, verbose.fig_output_file, ...
    fullfile(testFolder, 'report.pdf'));
verify_report_exists(testCase, verbose.fig_output_file);
end

function test_converts_legacy_postscript_to_pdf(testCase)
[verbose, testFolder] = create_test_figures(testCase, 'report.ps');

verifyWarning(testCase, ...
    @() tapas_physio_print_figs_to_file(verbose, testFolder), ...
    'tapas:physio:PostScriptNotSupported');
verbose = tapas_physio_print_figs_to_file_suppress_warning( ...
    verbose, testFolder);

verifyEqual(testCase, verbose.fig_output_file, ...
    fullfile(testFolder, 'report.pdf'));
verify_report_exists(testCase, verbose.fig_output_file);
end

function [verbose, testFolder] = create_test_figures(testCase, fileName)
testFolder = tempname;
mkdir(testFolder);
testCase.addTeardown(@() rmdir(testFolder, 's'));
figures = [figure('Visible', 'off'), figure('Visible', 'off')];
testCase.addTeardown(@() close(figures));
plot(axes('Parent', figures(1)), 1:3);
plot(axes('Parent', figures(2)), 3:-1:1);
verbose.level = 0;
verbose.fig_handles = figures;
verbose.fig_output_file = fileName;
end

function verbose = tapas_physio_print_figs_to_file_suppress_warning( ...
        verbose, testFolder)
warningState = warning('off', 'tapas:physio:PostScriptNotSupported');
cleanup = onCleanup(@() warning(warningState));
verbose = tapas_physio_print_figs_to_file(verbose, testFolder);
end

function verify_report_exists(testCase, fileName)
fileInfo = dir(fileName);
verifyTrue(testCase, isfile(fileName));
verifyGreaterThan(testCase, fileInfo.bytes, 0);
end
