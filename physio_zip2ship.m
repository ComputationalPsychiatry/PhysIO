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

zipFex = fullfile(currD, sprintf('PhysIOToolbox_r%s_examples.zip', srev));
zipFco = fullfile(currD, sprintf('PhysIOToolbox_r%s_code.zip', srev));
zip( zipFex, {fullfile(currD, 'examples'), fullfile(currD, 'README.txt')});
zip( zipFco, {fullfile(currD, 'code'), fullfile(currD, 'manual'),fullfile(currD, 'README.txt')});
