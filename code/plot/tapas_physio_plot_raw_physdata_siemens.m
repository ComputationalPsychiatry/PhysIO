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
set(fh, 'Name', stringTitle);

stringLegend = {'physio trigger on', 'physio trigger off', ...
    'assumed last sample of last scan volume', 'mean of channels'};

%% Plot each channel as subplot
nChannels = size(recordingChannels, 2);
layout = tiledlayout(fh, nChannels+1, 1, ...
    'TileSpacing', 'tight', 'Padding', 'compact');
if nargin >= 3
    % Reserve a footer outside the axes so timing stays visible during zoom.
    layout.OuterPosition = [0 0.10 1 0.90];
end
hs = gobjects(nChannels+1, 1);
for iChannel = 1:nChannels
    hs(iChannel) = nexttile(layout);
    channelData = recordingChannels(:,iChannel);
    channelOffset = 0;
    % Show small variations around large baselines without long tick labels.
    if abs(channelData(1)) > 100*max(1, max(channelData)-min(channelData))
        channelOffset = channelData(1);
    end
    plot(hs(iChannel), t, channelData-channelOffset);
    if channelOffset ~= 0
        hs(iChannel).YAxis.Exponent = 0;
        text(hs(iChannel), 0.99, 0.95, ...
            sprintf('offset: %+.0f', channelOffset), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', 'FontSize', 8, ...
            'BackgroundColor', 'white', 'Margin', 1);
    end
    ylabel(hs(iChannel), sprintf('channel %d', iChannel));
    hs(iChannel).XTickLabel = [];
end

%% Plot mean of all channels plus trigger stems in lowest subplot
hs(nChannels+1) = nexttile(layout);
axMean = hs(nChannels+1);
stem(axMean, cpulse_on, ampl*ones(size(cpulse_on)), 'g');
hold(axMean, 'on');
stem(axMean, cpulse_off, ampl*ones(size(cpulse_off)), 'r');
stem(axMean, t(stopSample), ampl , 'm');
plot(axMean, t, meanChannel);

ylabel(axMean, 'mean');

if ~isempty(recording_on)
    stem(axMean, recording_on, ampl*ones(size(recording_on)), 'k');
    stringLegend{end+1} = 'physio recording on';
end

if ~isempty(recording_off)
    stem(axMean, recording_off, ampl*ones(size(recording_off)), 'k');
    stringLegend{end+1} = 'physio recording off';
end

linkaxes(hs, 'x');
set(hs, 'FontSize', 11);
ytickformat(hs, '%.2f');
lg = legend(axMean, stringLegend, 'FontSize', 8, 'NumColumns', 2);
lg.Layout.Tile = 'north';
title(layout, stringTitle);
xlabel(axMean, 't (seconds)');

% Add text with physio logfile timing info
if nargin >=3
    textString{1} = sprintf('MDH/Scanner clock: start %.3f s | stop %.3f s | duration %.3f s', ...
        logFooter.StartTimeSecondsScannerClock, ...
        logFooter.StopTimeSecondsScannerClock, ...
        logFooter.StopTimeSecondsScannerClock - logFooter.StartTimeSecondsScannerClock);
    textString{2} = sprintf('MPCU/Physio Monitoring Unit clock: start %.3f s | stop %.3f s | duration %.3f s', ...
        logFooter.StartTimeSecondsRecordingClock, ...
        logFooter.StopTimeSecondsRecordingClock, ...
        logFooter.StopTimeSecondsRecordingClock - logFooter.StartTimeSecondsRecordingClock);
    annotation(fh, 'textbox', [0.02 0.005 0.96 0.085], ...
        'String', textString, 'FontSize', 8, 'Interpreter', 'none', ...
        'EdgeColor', 'none', 'VerticalAlignment', 'middle');
end
