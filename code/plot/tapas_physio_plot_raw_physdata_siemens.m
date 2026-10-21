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
% Reserve a shared left margin for the fixed-position channel labels.
layout.OuterPosition = [0.06 0 0.94 1];
if nargin >= 3
    % Reserve a footer outside the axes so timing stays visible during zoom.
    layout.OuterPosition = [0.06 0.10 0.94 0.90];
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
    hs(iChannel).YAxis.Exponent = 0;
    text(hs(iChannel), 0.99, 0.95, ...
        sprintf('offset: %+.12g', channelOffset), ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', 'FontSize', 8, ...
        'BackgroundColor', 'white', 'Margin', 1);
    ylabel(hs(iChannel), sprintf('channel %d', iChannel));
    hs(iChannel).XTickLabel = [];
end

%% Plot mean of all channels plus trigger stems in lowest subplot
hs(nChannels+1) = nexttile(layout);
axMean = hs(nChannels+1);
% Remove the baseline for display only; keep event stems on the trace scale.
meanOffset = meanChannel(1);
meanDisplay = meanChannel-meanOffset;
stemBase = min(meanDisplay);
stemHeight = max(meanDisplay) + 0.05*max(1, max(meanDisplay)-stemBase);
stem(axMean, cpulse_on, stemHeight*ones(size(cpulse_on)), 'g', ...
    'BaseValue', stemBase);
hold(axMean, 'on');
stem(axMean, cpulse_off, stemHeight*ones(size(cpulse_off)), 'r', ...
    'BaseValue', stemBase);
stem(axMean, t(stopSample), stemHeight, 'm', 'BaseValue', stemBase);
plot(axMean, t, meanDisplay);
axMean.YAxis.Exponent = 0;
text(axMean, 0.99, 0.95, sprintf('offset: %+.12g', meanOffset), ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', 'FontSize', 8, ...
    'BackgroundColor', 'white', 'Margin', 1);

ylabel(axMean, 'mean');

if ~isempty(recording_on)
    stem(axMean, recording_on, stemHeight*ones(size(recording_on)), 'k', ...
        'BaseValue', stemBase);
    stringLegend{end+1} = 'physio recording on';
end

if ~isempty(recording_off)
    stem(axMean, recording_off, stemHeight*ones(size(recording_off)), 'k', ...
        'BaseValue', stemBase);
    stringLegend{end+1} = 'physio recording off';
end

linkaxes(hs, 'x');
set(hs, 'FontSize', 11);
% Equal-width tick fields and fixed label positions keep all panels aligned.
ytickformat(hs, '%8.2f');
for iAxis = 1:numel(hs)
    hs(iAxis).YLabel.FontName = 'Helvetica';
    hs(iAxis).YLabel.FontSize = 10;
    hs(iAxis).YLabel.Units = 'normalized';
    hs(iAxis).YLabel.Position = [-0.10 0.5 0];
end
lg = legend(axMean, stringLegend, 'FontSize', 8, 'NumColumns', 2);
lg.Layout.Tile = 'north';
title(layout, stringTitle);
xlabel(axMean, 't (seconds)');

% Add text with physio logfile timing info
if nargin >=3
    formatTime = @(value) sprintf('%s (%.3f s)', ...
        char(duration(0, 0, value, 'Format', 'hh:mm:ss')), value);
    textString{1} = sprintf('MDH/Scanner clock: start %s | stop %s | duration %s', ...
        formatTime(logFooter.StartTimeSecondsScannerClock), ...
        formatTime(logFooter.StopTimeSecondsScannerClock), ...
        formatTime(logFooter.StopTimeSecondsScannerClock - logFooter.StartTimeSecondsScannerClock));
    textString{2} = sprintf('MPCU/Physio Monitoring Unit clock: start %s | stop %s | duration %s', ...
        formatTime(logFooter.StartTimeSecondsRecordingClock), ...
        formatTime(logFooter.StopTimeSecondsRecordingClock), ...
        formatTime(logFooter.StopTimeSecondsRecordingClock - logFooter.StartTimeSecondsRecordingClock));
    annotation(fh, 'textbox', [0.02 0.005 0.96 0.085], ...
        'String', textString, 'FontSize', 8, 'Interpreter', 'none', ...
        'EdgeColor', 'none', 'VerticalAlignment', 'middle');
end
