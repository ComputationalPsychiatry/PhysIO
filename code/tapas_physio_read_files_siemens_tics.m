function C = tapas_physio_read_files_siemens_tics(fileName, fileType)
% Reads _PULS, _RESP, _ECG, _Info-files from Siemens tics format with
% multiple numbers of columns and different column headers
%
%   output = tapas_physio_read_files_siemens_tics(input)
%
% IN
%   fileName    *.log from Siemens VD/VE tics file format
%   fileType    'ECG', 'PULS', 'RESP', 'Info'
%               If not specified, this is read from the last part of the
%               filename after the last underscore, e.g.
%               Physio_*_ECG.log -> log
% OUT
%
% EXAMPLE
%   tapas_physio_read_files_siemens_tics
%
%   See also
%
% Author: Lars Kasper
% Created: 2017-11-16
% Copyright (C) 2017 TNU, Institute for Biomedical Engineering,
%                    University of Zurich and ETH Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under the terms of the GNU General Public
% License (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id: teditRETRO.m 775 2015-07-17 10:52:58Z kasperla $

if nargin < 2
    % extract what is between last underscore and extension .log in
    % filename:
   fileType = regexp(fileName, '.*_([^_]*)\.log', 'tokens');
   if isempty(fileType)
       fileType = 'ECG';
   else
       fileType = upper(fileType{1});
   end

end