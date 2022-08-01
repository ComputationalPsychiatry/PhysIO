function [data_table, log_parts] = tapas_physio_siemens_line2table(lineData, cardiacModality)
% transforms data line of Siemens log file (before linebreak after 5003 
% signal for recording end) into table (sorting amplitude and trigger signals)
%
%   data_table = tapas_physio_siemens_line2table(input)
%
% IN
%   lineData        all recording amplitude samples are saved in first line
%                   of file. See also tapas_physio_read_physlogfiles_siemens_raw
%
% OUT
%   data_table      [nSamples,3] table of channel_1, channels_AVF and trigger
%                   signal with trigger codes:
%                   5000 = cardiac pulse on
%                   6000 = cardiac pulse off
%                   6002 = phys recording on
%                   6003 = phys recording off
%   log_parts       part of logfile according to markers by Siemens
%
% EXAMPLE
%   tapas_physio_siemens_line2table
%
%   See also

% Author: Lars Kasper
% Created: 2016-02-29
% Copyright (C) 2016 TNU, Institute for Biomedical Engineering,
%                    University of Zurich and ETH Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under the terms of the GNU General Public
% License (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.


iStartInfoRegion = strfind(lineData, '5002');
iEndInfoRegion = strfind(lineData, '6002');
iStartFooter = strfind(lineData, '5003');

% header marker identifier (data collection, signal source...)
logHeader = lineData(1:iStartInfoRegion(1)-1);

nInfoRegions = numel(iStartInfoRegion);

for iInfoRegion = 1:nInfoRegions
    logInfoRegions{iInfoRegion} = lineData((iStartInfoRegion(iInfoRegion)+3):(iEndInfoRegion(iInfoRegion)-1));
end


% Determine logfile version (which alters no of channels etc.) from header
%5002 LOGVERSION   1 6002%
%5002 LOGVERSION_RESP   3 6002%
%5002 LOGVERSION_PULS   3 6002%
logVersionCell = regexp(lineData, '5002 LOGVERSION[^0-9]*(\d+) 6002', 'tokens');
logVersion = str2double(logVersionCell{1}{1});


log_parts.logHeader = logHeader;
log_parts.logInfoRegions = logInfoRegion;
log_parts.logVersion = logVersion;

% last 6002 signals start of data logging
if ~isempty(iEndInfoRegion)
    iTrigger = iEndInfoRegion(end);
    % crop string after trigger
    lineData = lineData((iTrigger(end)+6):end);
    doCropLater = false;
else
    % crop first 4 values as in UPenn manual after conversion
    doCropLater = true;
end

data = textscan(lineData, '%d', 'Delimiter', ' ', 'MultipleDelimsAsOne',1);

if doCropLater
    % crop first 4 values;
    data{1} = data{1}(5:end);
end

% Remove the systems own evaluation of triggers.
cpulse  = find(data{1} == 5000);  % System uses identifier 5000 as trigger ON
cpulse_off = find(data{1} == 6000); % System uses identifier 5000 as trigger OFF
recording_on = find(data{1} == 6002);% Scanner trigger to Stim PC?
recording_off = find(data{1} == 5003);


% Filter the trigger markers from the ECG data
% Note: depending on when the scan ends, the last size(t_off)~=size(t_on).
iNonEcgSignals = [cpulse; cpulse_off; recording_on; recording_off];
codeNonEcgSignals = [5000*ones(size(cpulse)); ...
    6000*ones(size(cpulse_off)); ...
    6002*ones(size(recording_on))
    5003*ones(size(recording_off))];

% data_stream contains only the 2 ECG-channel time courses (with
% interleaved samples
data_stream = data{1};
data_stream(iNonEcgSignals) = [];

%iDataStream contains the indices of all true ECG signals in the full
%data{1}-Array that contains also non-ECG-signals
iDataStream = 1:numel(data{1});
iDataStream(iNonEcgSignals) = [];

nSamples = numel(data_stream);

switch upper(cardiacModality) % ecg has two channels, resp and puls only one
    case 'ECG'
        switch logVersion
            case 1
                nChannels = 2;
            case 3
                nChannels = 4;
        end
        nRows = ceil(nSamples/nChannels);
        
        % create a table with channel_1, channels_AVF and trigger signal in
        % different Columns
        % - iData_table is a helper table that maps the original indices of the
        % ECG signals in data{1} onto their new positions
        data_table = zeros(nRows,nChannels+1);
        iData_table = zeros(nRows,nChannels+1);
        
        for iChannel = 1:nChannels
            data_table(1:nRows,iChannel) = data_stream(iChannel:nChannels:end);
            iData_table(1:nRows,iChannel) = iDataStream(iChannel:nChannels:end);
        end
        
        % TODO: deal with mod(nSamples, nChannels) > 0 (incomplete data?)
        
        % now fill up nChannel+1. column with trigger data
        % - for each trigger index in data{1}, check where ECG data with closest
        % smaller index ended up in the data_table ... and put trigger code in
        % same row of that table
        nTriggers = numel(iNonEcgSignals);
        
        for iTrigger = 1:nTriggers
            % find index before trigger event in data stream and
            % detect it in table, look in last columns first, then go
            % backwards
            iRow = [];
            iChannel = nChannels;
            while isempty(iRow)
                
                iRow = find(iData_table(:,iChannel) == iNonEcgSignals(iTrigger)-1);
                iChannel = iChannel - 1;
            end
            
            data_table(iRow,nChannels+1) = codeNonEcgSignals(iTrigger);
        end
        
    case {'RESP', 'PPU'} % only one channel available, fill second row with zeros
        nRows = nSamples;
        
        % create a table with channel_1 and trigger signal in
        % different Columns
        % - iData_table is a helper table that maps the original indices of the
        % ECG signals in data{1} onto their new positions
        data_table = zeros(nRows,3);
        iData_table = zeros(nRows,3);
        
        data_table(1:nRows,1) = data_stream;
        iData_table(1:nRows,1) = iDataStream;
        
        % now fill up 3rd column with trigger data
        % - for each trigger index in data{1}, check where ECG data with closest
        % smaller index ended up in the data_table ... and put trigger code in
        % same row of that table
        nTriggers = numel(iNonEcgSignals);
        
        for iTrigger = 1:nTriggers
            % find index before trigger event in data stream and
            % detect it in table
            iRow = find(iData_table(:,1) == iNonEcgSignals(iTrigger)-1);
            if ~isempty(iRow)
                data_table(iRow,3) = codeNonEcgSignals(iTrigger);
            end
        end
    otherwise
        error('unknown cardiac/respiratory logging modality: %s', cardiacModality);
end