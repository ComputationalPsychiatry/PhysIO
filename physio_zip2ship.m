function physio_zip2ship(rev)
% deletes all unnecessary output files and zips code and examples folder
% separately 
%
%   output = physio_zip2ship(rev)
%
% IN
%   rev - revision number
% 
% OUT
%   PhysIOToolbox_r<rev>_code.zip
%   PhysIOToolbox_r<rev>_examples.zip
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
physio_cleanup_example_output_files();
currD = fileparts(mfilename('fullpath'));

exportD = fullfile(currD, sprintf('PhysIO_r%s', srev));
mkdir(exportD);
mkdir(exportD, 'manual');
copyfile(fullfile(currD, 'examples'), fullfile(exportD, 'examples'));
copyfile(fullfile(currD, 'code'), fullfile(exportD, 'code'));
copyfile(fullfile(currD, 'manual/*.pdf'), fullfile(exportD, 'manual'));
copyfile(fullfile(currD, 'README.txt'), exportD);

zipFex = fullfile(exportD, sprintf('PhysIOToolbox_r%s_examples.zip', srev));
zipFco = fullfile(exportD, sprintf('PhysIOToolbox_r%s_code.zip', srev));
zip( zipFex, {fullfile(exportD, 'examples'), fullfile(exportD, 'README.txt')});
zip( zipFco, {fullfile(exportD, 'code'), fullfile(exportD, 'manual'),fullfile(exportD, 'README.txt')});
