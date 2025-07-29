function pathToExamples = tapas_physio_download_example_data()
% Downloads example data for current version of PhysIO
%
%   pathToExamples = tapas_physio_download_example_data()
%
% IN
%
% OUT
%
% EXAMPLE
%   tapas_physio_download_example_data
%
%   See also
 
% Author:   Lars Kasper
% Created:  2025-06-17
% Copyright (C) 2025 TNU, Institute for Biomedical Engineering,
%                    University of Zurich and ETH Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under
% the terms of the GNU General Public License (GPL), version 3. You can
% redistribute it and/or modify it under the terms of the GPL (either
% version 3 or, at your option, any later version). For further details,
% see the file COPYING or <http://www.gnu.org/licenses/>.
 
currentRelease = tapas_physio_version();

semVersion = regexprep(currentRelease, '.*v', 'v');

% download current version of PhysIO examples, corresponding to code
% release version to temporary directory
urlZenodo = sprintf('https://zenodo.org/records/15579087/files/ComputationalPsychiatry/PhysIO-Examples-%s.zip', semVersion);
tempZipFilePath = [tempname '.zip'];  % tempname is matlab inbuilt
fprintf('Downloading Example Data for PhysIO version %s into PhysIO/examples folder...\n', semVersion);
fprintf('This may take a few minutes (50 MB)\n')
websave(tempZipFilePath, urlZenodo);
fprintf('Done. \n\n')

% unzip examples into canonical PhysIO Examples folder
doVerifyPath = false;
pathToExamples = tapas_physio_get_path_examples([], ...
    doVerifyPath);

unzip(tempZipFilePath, pathToExamples);

% Cleanup: examples are in a subfolder
% ComputationalPsychiatry-PhysIO-Examples-<GitCommitHash>, move them
% directly into examples
dirPhysioExamples = dir([pathToExamples '/ComputationalPsychiatry-PhysIO-Examples-*']).name;
movefile(fullfile(pathToExamples, dirPhysioExamples, "*"), pathToExamples)
rmdir(fullfile(pathToExamples, dirPhysioExamples), 's')
