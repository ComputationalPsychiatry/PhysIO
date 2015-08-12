function physio_zip2ship(rev, zipWhat)
% deletes all unnecessary output files and zips code and examples folder
% separately
%
%   output = physio_zip2ship(rev, zipWhat)
%
% IN
%   rev     revision number as a string or integer value
%   zipWhat  determines which part of the code will be zipped
%           'withSiemensTics'       copy proprietary specs of Siemens Tics
%                                   logfiles as well
%           'withoutSiemensTics'    (default)
%                                   do not copy Siemens tics code
%           'all'                   one big file with code and examples
% OUT
%   PhysIOToolbox_r<rev>_code.zip
%   PhysIOToolbox_r<rev>_examples.zip
%   PhysIOToolbox_r<rev>_all.zip    for zipWhat 'all', code and
%                                   examples are combined in this logfile
%
% EXAMPLE
%   physio_zip2ship(193)
%
%   See also physio_cleanup_example_output_files
%
% Author: Lars Kasper
% Created: 2013-05-10
% Copyright (C) 2013 TNU, Institute for Biomedical Engineering, University of Zurich and ETH Zurich.
%
% This file is part of the TNU CheckPhysRETROICOR toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id$

if ~nargin
    srev = datestr(now, 'yyyy_mm_dd_HHMMSS'); %TODO: somehow parse physio_new to get a better way to detect revision
else
    srev = sprintf('%d', rev');
end

if nargin < 2
    zipWhat = 'all';
end

% delete created output files, e.g. physio.mat, multiple_regressors.txt
physio_cleanup_example_output_files();

% update job.mat and matlab_script from job.m, remove absolute paths
physio_update_examples();
currD = fileparts(mfilename('fullpath'));

exportD = fullfile(currD, sprintf('PhysIO_r%s', srev));
mkdir(exportD);

% copy and zip examples with manual
copyfile(fullfile(currD, 'README.txt'), exportD);

dirExamples = {
    'GE/PPU3T'
    'Philips/ECG3T'
    'Philips/ECG7T'
  %  'Philips/PPU3T'
    'Siemens/ECG3T'
    };

nExamples = numel(dirExamples);
for iExample = 1:nExamples
    mkdir(fullfile(exportD, 'examples', dirExamples{iExample}));
    copyfile(fullfile(currD, 'examples', dirExamples{iExample}), ...
        fullfile(exportD, 'examples',dirExamples{iExample}));
end
mkdir(exportD, 'manual');
copyfile(fullfile(currD, 'manual/*.pdf'), fullfile(exportD, 'manual'));

zipF = fullfile(currD, sprintf('PhysIOToolbox_r%s_examples.zip', srev));
zip(zipF, exportD);

copyfile(fullfile(currD, 'code'), fullfile(exportD, 'code'));

% copy and zip all
switch zipWhat
    case 'all'
        zipF = fullfile(currD, sprintf('PhysIOToolbox_r%s_all.zip', srev));
        zip(zipF, exportD);
    case 'withoutSiemensTics'
        fileArraySiemensTics = {
            'tapas_physio_create_scan_timing_from_tics_siemens.m'
            'tapas_physio_read_physlogfiles_siemens_tics.m'
            };
        nFiles = numel(fileArraySiemensTics);
        for iFile = 1:nFiles
            delete(fullfile(exportD, 'code', fileArraySiemensTics{iFile}));
        end
end

% delete large data files and save examples wo data
threshSize = 1e6;

for iExample = 1:nExamples
    currDir = fullfile(exportD, 'examples', dirExamples{iExample});
    fileList=dir(currDir);
    nFiles = numel(fileList);
    for iFile = 1:nFiles
        if fileList(iFile).bytes > threshSize
            delete(fullfile(currDir, fileList(iFile).name));
            
            % create empty file with download note
            fileDummy = fullfile(currDir, [fileList(iFile).name ...
                '_DOWNLOAD_FROM_TNU-ZURICH.ORG_SOFTWARE_TAPAS-DATA']);
            fclose(fopen(fileDummy, 'w'));
            
        end
    end
end

zipF = fullfile(currD, sprintf('PhysIOToolbox_r%s_code_examples_wo_data.zip', srev));
zip(zipF, exportD);


% copy and zip code and manual
rmdir(fullfile(exportD, 'examples'),'s');
rmdir(fullfile(exportD, 'manual'),'s')
mkdir(exportD, 'manual');
copyfile(fullfile(currD, 'manual/QuickStart*.pdf'), fullfile(exportD, 'manual'));

zipF = fullfile(currD, sprintf('PhysIOToolbox_r%s_code.zip', srev));
zip(zipF, exportD);

