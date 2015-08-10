function output = physio_update_examples(input)
% Copies current tapas_physio_examples (_spm_job.m) to example folder, remove absolute
% paths, create _spm_job.mat and _matlab_script from it
%
%   output = physio_update_examples(input)
%
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

pathRoot = fileparts(mfilename('fullpath'));
pathExamples = fullfile(oathRoot, 'examples');
pathCode = fullfile(pathRoot, 'code');

pfxExample = 'tapas_physio_example_';

dirExamples = {
    'Philips/ECG3T'
    'Philips/ECT7T'
    'Philips/PPU3T'
    'Siemens/ECG3T'
    'GE/PPU3T'
    };


% copy file to example-folder, remove(?) tapas_physio

% remove absolute folder definitions, also for SPM (?)

% save accompanying .mat-spm job and matlab-script
physio = tapas_physio_save_batch_mat_script(fileBatchM);
