function isPhysioCorrectlyInitialized = tapas_physio_init()
% Initializes TAPAS by checking that all folders are at the right
%
%    isPhysioCorrectlyInitialized = tapas_physio_init()
%
% IN
%
% OUT
%
% EXAMPLE
%   tapas_physio_init()
%
%   See also
%
% Author: Lars Kasper
% Created: 2018-02-17
% Copyright (C) 2018 TNU, Institute for Biomedical Engineering,
%                    University of Zurich and ETH Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under the terms of the GNU General Public
% License (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
disp(' ');
disp(' ');
disp('====================================================================')
disp(' ');
disp('          This is the TAPAS PhysIO Toolbox Version R2017.3');
disp(' ');
disp('====================================================================')
disp(' ');
disp('Checking Matlab and PhysIO paths now');

[isPhysioVisibleForSpmBatchEditor, pathSpm, pathPhysIO] = ...
    tapas_physio_check_spm_batch_editor_integration();

isPhysioCorrectlyInitialized = ~isempty(pathSpm) && ~isempty(pathPhysIO) && ...
    isPhysioVisibleForSpmBatchEditor;

if isPhysioCorrectlyInitialized
    disp('Success: PhysIO successfully installed, integration with Batch Editor possible.')
    fprintf('Updating SPM batch editor information...')
    spm_jobman('initcfg');
    fprintf('done.\n\n');
else
    if isempty(pathPhysIO)
        pathPhysIOHere = fileparts(mfilename('fullpath'));
        warning(['\n PhysIO not setup correctly, add %s to the path, e.g., via ' ...
            '\n         addpath %s'],  pathPhysIOHere)
    else
        warning(['\n PhysIO integration with SPM not setup. \n' ...
            ' PhysIO will run in Matlab, but not with SPM Batch Editor GUI. \n' ...
            ' Follow instructions above to link the paths for batch editor use.'])
    end
end