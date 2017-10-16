function [VOLLOCS, LOCS, verbose] = ...
    tapas_physio_create_scan_timing_from_tics_siemens(t, t_start, log_files, verbose)
% Creates locations of scan volume and slice events in physiological time series vector from Siemens Tics-file
% (<date>_<time>_AcquisitionInfo*.log)
%
% [VOLLOCS, LOCS, verbose] = ...
%    tapas_physio_create_scan_timing_from_tics_siemens(t, log_files, verbose)
%
% IN
%   t           - timing vector of physiological logfiles
%   t_start       offset to t indicating start of phys log file
%   log_files   - structure holding log-files, here:
%                 .scan_timing - <date>_<time>_AcquisitionInfo*.log
%                                Siemens Acquisition info logfile, e.g.
%                                 Volume_ID Slice_ID AcqTime_Tics
%                                 0 0 22953292
%                                 0 1 22954087
%                                 0 2 22953690
%                                 1 0 22954700
%                                 1 1 22955495
%                                 1 2 22955097
%                                 2 0 22956136
% OUT
%           VOLLOCS         - locations in time vector, when volume scan
%                             events started
%           LOCS            - locations in time vector, when slice or volume scan
%                             events started
% EXAMPLE
%   [VOLLOCS, LOCS] = tapas_physio_create_nominal_scan_timing(t, sqpar);
%
%   See also
%
% Author: Lars Kasper
% Created: 2014-09-08
% Copyright (C) 2014 Institute for Biomedical Engineering, ETH/Uni Zurich.
%
% This file is part of the PhysIO toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id$
DEBUG = verbose.level >=3;

fid = fopen(log_files.scan_timing);

% Determine number of header lines by searching for the column header line,
% which has both Volume and Slice as a keyword in it
haveFoundColumnHeader = false;
nHeaderLines = 0;
while ~haveFoundColumnHeader
    nHeaderLines = nHeaderLines + 1;
    strLine = fgets(fid);
    haveFoundColumnHeader = any(regexp(upper(strLine), 'VOLUME.* *SLICE'));
end
fclose(fid);

nColumns = numel(regexp(strLine, ' *')) + 1; % columns are separated by arbitrary number of spaces
nHeaderLines = nHeaderLines + 1; % since empty line after header

fid = fopen(log_files.scan_timing);
switch nColumns
    case 3
        C = textscan(fid, '%d %d %d', 'HeaderLines', 1); % no empty line here...
    case 4
        % 4 columns variant: VOLUME   SLICE   ACQ_START_TICS  ACQ_FINISH_TICS
        C           = textscan(fid, '%d %d %d %d', 'HeaderLines', nHeaderLines);
    case 5
        % NOTE: a Tic of 0 is extremely unlikely and indicates a missed header
        % line, so go for the 5-columns variant:
        % VOLUME   SLICE   ACQ_START_TICS  ACQ_FINISH_TICS  ECHO
        C           = textscan(fid, '%d %d %d %d %d', 'HeaderLines', nHeaderLines);
end


dtTicSeconds = 2.5e-3;

idVolumes       = C{1};
idSlices        = C{2};

% HACK: take care of interleaved acquisition; multiband?
ticsAcq         = sort(double(C{3}));
tAcqSeconds     = ticsAcq*dtTicSeconds - t_start; % relative timing to start of phys logfile

% find time in physiological log time closest to time stamp of acquisition
% for each time
if exist('knnsearch', 'file')
    LOCS = knnsearch(t, tAcqSeconds); % Matlab stats toolbox next neighbor search
else
    % slower linear search, without the need for knnsearch
    LOCS = zeros(size(idSlices));
    iSearchStart = 1;
    for iLoc = 1:numel(LOCS)
        if ~mod(iLoc,1000), fprintf('%d/%d\n',iLoc,numel(LOCS));end
        [~, iClosestTime] = min(abs(t(iSearchStart:end) - tAcqSeconds(iLoc)));
        
        % only works for ascendingly sorted ticsAcq!
        LOCS(iLoc) = iSearchStart - 1 + iClosestTime;
        iSearchStart = LOCS(iLoc);
    end
end

% extract start times of volume by detecting index change volume id
indVolStarts = [1; find(diff(idVolumes) > 0) + 1]; 
VOLLOCS = LOCS(indVolStarts);

if DEBUG
   verbose.fig_handles(end+1) = tapas_physio_get_default_fig_params();
   stringTitle = 'Extracted Sequence Timing Siemens';
   set(gcf, 'Name', stringTitle);
   stem(t(LOCS), ones(size(LOCS))); hold all;
   stem(t(VOLLOCS), 1.5*ones(size(VOLLOCS)));
   legend('slice start events', 'volume start events');
   xlabel('t (seconds)');
   title(stringTitle);
end