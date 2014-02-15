function fh = tapas_physio_plot_raw_physdata_diagnostics(cpulse, yResp, thresh_cardiac)
% plots diagnostics for raw physiological time series as monitoried by the
% MR scanner breathing belt/ECG
%
% Author: Lars Kasper
%
% Copyright (C) 2013, Institute for Biomedical Engineering, ETH/Uni Zurich.
%
% This file is part of the PhysIO toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id$

% cardiac analysis of heartbeat rates

hasCardiacData = ~isempty(cpulse);
hasRespData = ~isempty(yResp);
fh = tapas_physio_get_default_fig_params();
set(fh, 'Name','Diagnostics raw phys time series');
ah = subplot(2,1,1);

if hasCardiacData
    percentile = thresh_cardiac.percentile;
    upperThresh = thresh_cardiac.upperThresh;
    lowerThresh = thresh_cardiac.lowerThresh;
    [outliersHigh,outliersLow,fh] = tapas_physio_cardiac_detect_outliers(cpulse, percentile, upperThresh, lowerThresh, ah);
end
title( 'temporal lag between subsequent heartbeats (seconds)');

% histogram of breathing amplitudes
subplot(2,1,2);

if hasRespData
    nBins = min(length(unique(yResp)), floor(length(yResp)/100));
    hist(yResp, nBins);
end
title('histogram of breathing belt amplitudes');
end
