function physio_cleanup_example_output_files()
% deletes all saved figures and multiple_regressors and .asv-files in
% examples subdirectory to make it ready for zipping and shipping
%
%  physio_cleanup_example_output_files()
%
% IN
%
% OUT
%
% EXAMPLE
%   cleanup_example_output_files.m
%
%   See also physio_zip2ship.m
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

currD = fileparts(mfilename('fullpath'));
pfxD = {'examples/GE';'examples/Philips'; 'code'; 'manual'};

for iPfx = 1:length(pfxD)
    ds = dir(fullfile(currD, pfxD{iPfx}));
    ds = ds(3:end);
    
    for n = 1:length(ds);
        d = fullfile(currD,pfxD{iPfx},ds(n).name);
        disp(['Deleting output files in ' d]);
        delete(fullfile(d, 'multiple_regressors.*'));
        delete(fullfile(d, 'PhysIO_output*'));
        delete(fullfile(d, '*.asv'));
        %     delete(fullfile(d, '*.ps'));
        %     delete(fullfile(d, '*.tif*'));
        %     delete(fullfile(d, '*.jpg'));
        %     delete(fullfile(d, '*.jpeg'));
        %     delete(fullfile(d, '*.fig'));
        %     delete(fullfile(d, '*.jpg'));
    end
end
