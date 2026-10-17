function fh = tapas_physio_plot_raw_physdata_siemens(dataPhysio, stringModality, logFooter)
% plots cardiac data as extracted from Siemens log file
%
%   fh = tapas_physio_plot_raw_physdata_siemens(dataCardiac)
%
% IN
%   dataPhysio      output struct from tapas_physio_siemens_table2cardiac
%   stringModality  for title display e.g., log file extension(ECG/PULS/RESP), 
%   logFooter       for timing display of MDH/MPCU time
% OUT
%
% EXAMPLE
%   tapas_physio_plot_raw_physdata_siemens
%
%   See also tapas_physio_read_physlogfiles_siemens tapas_physio_siemens_table2cardiac

% Author: Lars Kasper
% Created: 2016-02-29
% Copyright (C) 2016 TNU, Institute for Biomedical Engineering,
%                    University of Zurich and ETH Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under the terms of the GNU General Public
% License (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.


tapas_physio_strip_fields(dataPhysio);

stringTitle = 'Read-In: Raw Siemens physlog data';

if nargin >=2
    stringTitle = [stringTitle ' (' stringModality ')'];
end

fh = tapas_physio_get_default_fig_params();
set(gcf, 'Name', stringTitle);

stringLegend = {'raw channel data', 'physio trigger on', 'physio trigger off', ...
    'assumed last sample of last scan volume'};

nChannels = size(recordingChannels, 2);
for iChannel = 1:nChannels
    hs(iChannel) = subplot(nChannels+1, 1, iChannel);
    plot(t, recordingChannels(:,iChannel));
    title(sprintf('channel %d', iChannel));
end

hs(nChannels+1) = subplot(nChannels+1, 1, nChannels+1);
stem(cpulse_on, ampl*ones(size(cpulse_on)), 'g'); hold all;
stem(cpulse_off, ampl*ones(size(cpulse_off)), 'r');
stem(t(stopSample), ampl , 'm');
plot(t, meanChannel);

title('mean of channels');

if ~isempty(recording_on)
    stem(recording_on, ampl*ones(size(recording_on)), 'k');
    stringLegend{end+1} = 'physio recording on';
end

if ~isempty(recording_off)
    stem(recording_off, ampl*ones(size(recording_off)), 'k');
    stringLegend{end+1} = 'physio recording off';
end

linkaxes(hs, 'x')
legend(stringLegend);
tapas_physio_suptitle(stringTitle);
xlabel('t (seconds)');

% Add text with physio logfile timing info
if nargin >=3
    textString{1} = sprintf('StartTime (MDH/Scanner clock): %f s', logFooter.StartTimeSecondsScannerClock);
    textString{2} = sprintf('StopTime (MDH/Scanner clock): %f s (delta: %f s)', ...
        logFooter.StopTimeSecondsScannerClock, ...
        logFooter.StopTimeSecondsScannerClock - logFooter.StartTimeSecondsScannerClock);
    textString{3} = sprintf('StartTime (MPCU/Physio Monitoring Unit clock): %f s', logFooter.StartTimeSecondsRecordingClock);
    textString{4} = sprintf('StopTime (MPCU/Physio Monitoring Unit clock): %f s (delta: %f s)', ...
        logFooter.StopTimeSecondsRecordingClock, ...
        logFooter.StopTimeSecondsRecordingClock - logFooter.StartTimeSecondsRecordingClock);
    text(0,0, textString);
end