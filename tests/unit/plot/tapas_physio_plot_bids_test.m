function tests = tapas_physio_plot_bids_test()
% Checks BIDS read-in data and events passed to existing review plots
%
% Compares the PPU example with its JSON/TSV source data and checks that
% the plotted cardiac trace uses the same samples and time coordinates.
% A separate fixture checks detected peaks stored in cardiac_pulse.
%
% EXAMPLE
%   results = runtests('tapas_physio_plot_bids_test');
%
%   See also tapas_physio_plot_bids tapas_physio_readin_bids_test

tests = functiontests(localfunctions);
end

%% Locate existing example data, needed for the PPU test
function setupOnce(testCase)
% Use version-matched Zenodo paths and verify local data without downloading.
doUseZenodoPaths = true; doVerifyPath = true; doDownloadData = false;
[testCase.TestData.pathExamples, ...
    testCase.TestData.pathTestReferenceResults] = ...
    tapas_physio_get_paths_for_tests(...
    doUseZenodoPaths, doVerifyPath, doDownloadData);
end

%% Compare read-in and plotted PPU traces with the BIDS source files
function test_ppu_matches_source(testCase)
% Pass the JSON filename to test inference of the compressed TSV.
% Create hidden figures and close only those created by this test.
filename = fullfile(testCase.TestData.pathExamples, 'BIDS', 'PPU3T', ...
    'sub-s999_task-random_run-99_physio.json');
defaults = tapas_physio_new();
defaults.verbose.show_figs = false;
[verbose, actual, sqpar] = tapas_physio_plot_bids(...
    filename, [], defaults.verbose);
testCase.addTeardown(@() close(verbose.fig_handles));
% Read the TSV independently of the PhysIO reader. Unzip into a
% temporary folder so that the example data remain unchanged.
folder = tempname;
mkdir(folder);
testCase.addTeardown(@() rmdir(folder, 's'));
files = gunzip(strrep(filename, '.json', '.tsv.gz'), folder);
expected = readmatrix(files{1}, 'FileType', 'text');
metadata = jsondecode(fileread(filename));
% BIDS sample times are relative to the first volume, including
% the negative StartTime for samples recorded before acquisition.
expectedTime = metadata.StartTime + ...
    (0:size(expected,1)-1)' / metadata.SamplingFrequency;
testCase.verifyEqual(actual.t, expectedTime, 'AbsTol', 1e-12);
testCase.verifyEqual(actual.c, expected(:,1));
testCase.verifyEqual(actual.r, expected(:,2));

% Without additional acquisition parameters, use one slice per
% volume and unchanged respiration for the display fields.
testCase.verifyEqual(actual.spulse, actual.svolpulse);
testCase.verifyEqual(actual.fr, actual.r);
testCase.verifyEqual(sqpar.Nslices, 1);
testCase.verifyNumElements(verbose.fig_handles, 2);

% Check the plotted coordinates as well as the returned data:
% select the solid red cardiac trace from the existing plot.
lines = findobj(verbose.fig_handles(1), 'Type', 'line', ...
    'Color', [1 0 0], 'LineStyle', '-');
testCase.assertNumElements(lines, 1);
testCase.verifyEqual(lines.YData(:), actual.c);
testCase.verifyEqual(lines.XData(:), actual.t);
end

%% Preserve exported cardiac peaks when no scan triggers are present
function test_derivative_retains_cardiac_peaks(testCase)
% Create a small BIDS derivative with a separate cardiac_pulse
% column, as written by tapas_physio_write2bids at step 4.
% The all-zero trigger column represents no recorded scan events.
folder = tempname;
mkdir(folder);
testCase.addTeardown(@() rmdir(folder, 's'));
filename = fullfile(folder, 'test_physio');
% Peaks at samples 20 and 60 correspond to 0.9 and 4.9 seconds
% with StartTime = -1 s and SamplingFrequency = 10 Hz.
data = [(1:100)', sin((1:100)'), zeros(100,1), ...
    double(ismember((1:100)', [20 60]))];
writematrix(data, [filename '.tsv'], 'FileType', 'text', ...
    'Delimiter', 'tab');
metadata = struct('StartTime', -1, 'SamplingFrequency', 10, ...
    'Columns', {{'cardiac','respiratory','trigger','cardiac_pulse'}});
fid = fopen([filename '.json'], 'w');
fprintf(fid, '%s', jsonencode(metadata));
fclose(fid);
defaults = tapas_physio_new();
defaults.verbose.show_figs = false;
[verbose, actual] = tapas_physio_plot_bids(...
    [filename '.tsv'], [], defaults.verbose);
testCase.addTeardown(@() close(verbose.fig_handles));
% Read the uncompressed TSV directly and verify that peak times
% survive while the absent volume events remain empty.
testCase.verifyEqual(actual.cpulse, [0.9; 4.9], 'AbsTol', 1e-12);
testCase.verifyEmpty(actual.svolpulse);
end
