function pathTestReferenceResults = tapas_physio_download_test_reference_results()
% Downloads reference results for tests of current version of PhysIO
%
%   pathTestReferenceResults = tapas_physio_download_test_reference_results()
%
% IN
%
% OUT
%
% EXAMPLE
%   tapas_physio_download_test_reference_results
%
%   See also
 
% Author:   Lars Kasper
% Created:  2025-07-29
% Copyright (C) 2025 TNU, Institute for Biomedical Engineering,
%                    University of Zurich and ETH Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under
% the terms of the GNU General Public License (GPL), version 3. You can
% redistribute it and/or modify it under the terms of the GPL (either
% version 3 or, at your option, any later version). For further details,
% see the file COPYING or <http://www.gnu.org/licenses/>.
 
currentRelease = tapas_physio_version();

% TODO: make this a new version / find latest existing one online
%semVersion = regexprep(currentRelease, '.*v', 'v');
semVersion = 'v9.0.3';

% download current version of PhysIO examples, corresponding to code
% release version to temporary directory
urlZenodo = sprintf('https://zenodo.org/records/16579519/files/ComputationalPsychiatry/PhysIO-Test-Reference-Results-%s.zip', semVersion);
tempZipFilePath = [tempname '.zip'];  % tempname is matlab inbuilt
fprintf('Downloading Test Reference Results for PhysIO version %s into PhysIO/test-reference-results folder...\n', semVersion);
fprintf('This may take a few minutes (250 MB)\n')
websave(tempZipFilePath, urlZenodo);
fprintf('Done. \n\n')

% unzip examples into canonical PhysIO Examples folder
doVerifyPath = false;

% make sure subfolders are in path before calling any other functions


pathTestReferenceResults = tapas_physio_get_path_test_reference_results([], ...
    doVerifyPath);

unzip(tempZipFilePath, pathTestReferenceResults);

% Cleanup: examples are in a subfolder
% ComputationalPsychiatry-PhysIO-Examples-<GitCommitHash>, move them
% directly into examples
dirPhysioExamples = dir([pathTestReferenceResults '/ComputationalPsychiatry-PhysIO-Test-Reference-Results-*']).name;
movefile(fullfile(pathTestReferenceResults, dirPhysioExamples, "*"), pathTestReferenceResults)
rmdir(fullfile(pathTestReferenceResults, dirPhysioExamples), 's')
