function dirExamples = physio_update_examples()
% Copies current tapas_physio_examples from code-folder (_spm_job.m) 
% to example folder, remove absolute paths, create _spm_job.mat and _matlab_script from it
%
%   dirExamples = physio_update_examples()
%
% NOTE: The replacement of the folder names was only tested for Mac/Unix
%       systems
% IN
%
% OUT
%
% EXAMPLE
%   physio_update_examples
%
%   See also
%
% Author: Lars Kasper
% Created: 2015-08-10
% Copyright (C) 2015 TNU, Institute for Biomedical Engineering,
%                    University of Zurich and ETH Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under the terms of the GNU General Public
% License (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id$

pathRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
pathExamples = fullfile(pathRoot, 'private', 'examples');
pathCode = fullfile(pathRoot, 'public', 'code');
pathSpm = fileparts(which('spm'));

if isempty(pathSpm)
    error('Add SPM to your path');
else
    %addpath(genpath(fullfile(pathSpm, 'matlabbatch')));
    spm_jobman('initcfg');
end

pfxExample = 'tapas_physio_example_';

dirExamples = {
    'Philips/ECG3T'
    'Philips/ECG3T_V2'
    'Philips/ECG7T'
    'GE/PPU3T'
    'Siemens/ECG3T'
 %   'Philips/PPU3T'
    };

nExamples = numel(dirExamples);

for iExample = 1:nExamples
    dirExample = dirExamples{iExample};
    pathExample = fullfile(pathExamples, dirExample);
    fileJobM = [regexprep(lower(dirExample), '/', '_') '_spm_job.m'];
    fileJobSource = fullfile(pathCode, [pfxExample, fileJobM]);
    fileJobTarget = fullfile(pathExample, fileJobM);
    
    % open new file to write in example-folder, remove tapas_physio from
    % name
    
    % write file new, remove absolute folder definitions, also for SPM (?)
    fidSource = fopen(fileJobSource, 'r');
    fidTarget = fopen(fileJobTarget, 'w+');
    while 1
        tline = fgetl(fidSource);
        if ischar(tline)
            % replace absolute path
            tline = regexprep(tline, [pathExample '[/]*'], '');
            fprintf(fidTarget, '%s\n', tline);
        else
            break
        end
    end
    fclose(fidTarget);
    fclose(fidSource);
    
    % save accompanying .mat-spm job and matlab-script
    physio = tapas_physio_save_batch_mat_script(fileJobTarget);
end