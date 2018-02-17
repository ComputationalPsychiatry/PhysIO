function cmdString = tapas_physio_create_spm_toolbox_link()
% Creates a symbolik link of PhysIO/code folder to subfolder SPM/toolbox/PhysIO
% to make toolbox visible to SPM Batch editor
%
%   tapas_physio_create_spm_toolbox_link()
%
% IN
%
% OUT
%   cmdString   string of executed command
%
% EXAMPLE
%   tapas_physio_create_spm_toolbox_link
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
cmdString = '';

pathPhysIO = fileparts(mfilename('fullpath'));

pathSpm = fileparts(which('spm'));

if isempty(pathSpm)
    warning('SPM folder not found. Could not create symbolink link to PhysIO Toolbox');
else
    pathLinkPhysIOSPM = fullfile(pathSpm, 'toolbox', 'PhysIO');
    if ~exist(pathLinkPhysIOSPM, 'dir')
        % Create Link, OS dependent
        if ispc
            % indeed the other way around than in Linux/Mac
            cmdString = sprintf('mklink /D %s %s', pathLinkPhysIOSPM, pathPhysIO);
        else %unix/Max
            cmdString = sprintf('ln -s %s %s', pathPhysIO, pathLinkPhysIOSPM);
        end
        system(cmdString);
    end
end