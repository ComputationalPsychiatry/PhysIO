function plot_raw_physdata_diagnostics(t, tCardiac, yResp)
% plots diagnostics for raw physiological time series as monitoried by the
% MR scanner breathing belt/ECG
%
% Author: Lars Kasper
%
% Copyright (C) 2013, Institute for Biomedical Engineering, ETH/Uni Zurich.
%
% This file is part of the TNU CheckPhysRETROICOR toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id$

fh = get_default_fig_params();
set(fh, 'Name','Diagnostics raw phys time series');
subplot(2,1,1);
dt = diff(tCardiac);

plot(tCardiac(2:end), dt);
xlabel('t (seconds)');
ylabel('lag between heartbeats (seconds)');
title('temporal lag between heartbeats');


percentile = 0.8;
deviationPercent = 60;

prctileValue = my_prctile(dt, percentile);

if max(dt) > (1+deviationPercent/100)*prctileValue
    text(t( find(dt==max(dt))+1 ),max(dt),...
        {'Warning: There seem to be skipped heartbeats R-waves in the scanner-log', ...
        'rerun read\_physlog with ECG\_min set to 1'}, ...
        'Color', [1 0 0])
end

subplot(2,1,2);
nBins = min(length(unique(yResp)), floor(length(yResp)/100));
hist(yResp, nBins);
title('histogram of breathing belt amplitudes');
end