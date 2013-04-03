function [c, r, t, cpulse] = read_physlogfiles(files, cardiac_modality)
% reads out physiological time series and timing vector depending on the
% MR scanner vendor and the modality of peripheral cardiac monitoring (ECG
% or pulse oximetry)
%
%   [cpulse, rpulse, t, c] = read_physlogfiles(logfile, vendor, cardiac_modality)
%
% IN
%   files   is a structure containing the following filenames (with full
%           path)
%       .vendor             'Philips', 'GE' or 'Siemens', depending on your
%                           MR Scanner system
%       .log_cardiac        contains ECG or pulse oximeter time course
%                           for Philips: 'SCANPHYSLOG<DATE&TIME>.log';
%                           can be found on scanner in G:/log/scanphyslog-
%                           directory, one file is created per scan, make sure to take
%                           the one with the time stamp corresponding to your PAR/REC
%                           files
%       .log_respiration    contains breathing belt amplitude time course
%                           for Philips: same as .log_cardiac
%   cardiac_modality    'ECG' for ECG, 'OXY' for pulse oximetry
%
% OUT
%   cpulse              time events of R-wave peak in cardiac time series (seconds)
%   r                   respiratory time series
%   t                   vector of time points (in seconds)
%   c                   cardiac time series (ECG or pulse oximetry)
%
% EXAMPLE
%   [ons_secs.cpulse, ons_secs.rpulse, ons_secs.t, ons_secs.c] =
%   read_physlogfiles(logfile, vendor, cardiac_modality);
%
%   See also main_create_retroicor_regressors
%
% Author: Lars Kasper
% Created: 2013-02-16
% Copyright (C) 2013, Institute for Biomedical Engineering, ETH/Uni Zurich.
%
% This file is part of the TNU CheckPhysRETROICOR toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id$

switch lower(files.vendor)
    case 'philips'
        % everything stored in 1 logfile
        if ~isfield(files, 'log_cardiac') || isempty(files.log_cardiac)
            logfile = files.log_respiration;
        else
            logfile = files.log_cardiac;
        end
        [c, r, t, cpulse] = read_physlogfiles_philips(logfile, cardiac_modality);
    case 'ge'
        [c, r, t, cpulse] = read_physlogfiles_GE(files);
    case 'siemens'
        disp('Ask the FIL about it...');
end
