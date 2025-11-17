function  fh = tapas_physio_plot_movement_outliers_fd(rp, quality_measures, ...
    censoring, censoring_threshold, verbose, TR)
% Plots movement outliers (for censoring), based on framewise displacement
% (FD), as computed by Power et al., 2012. Also the plotting style is based
% on Power et al., 2017
%
%  fh = tapas_physio_plot_movement_outliers_fd(rp, quality_measures, ...
%           censoring, censoring_threshold, TR)
%
% IN
%   quality_measures    - output of tapas_physio_get_movement_quality_measures
%   censoring           - output tapas_physio_create_movement_regressors,
%                         default: []
%                         censoring = struct('nOutliers', nOutliers, 'R_outlier', R_outlier, ...
%                           'iOutlierTrans', iOutlierTrans, 'iOutlierRot', iOutlierRot, ...
%                           'iOutlierArray', iOutlierArray);
%
%   censoring threshold - for horizontal line indicating censoring
%                         default: from tapas_physio_new()
%   verbose             - for saving figure handlles
%                         default: []
%   TR or t               - repetion time (in seconds) for fMRI volumes
%                           or correct time vector (nScans,1)
%                         default: [] (plotted in multiples of volume)
%
% OUT
%
% EXAMPLE
%   tapas_physio_plot_movement_outliers_fd
%
%   See also tapas_physio_get_movement_quality_measures tapas_physio_create_movement_regressors

% Author: Lars Kasper
% Created: 2018-02-21
% Copyright (C) 2018 TNU, Institute for Biomedical Engineering,
%                    University of Zurich and ETH Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under the terms of the GNU General Public
% License (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%

if nargin < 6
    TR = [];
end

hasTR = ~isempty(TR);

if nargin >= 5
    % If verbose is passed as argument (from updated tapas_physio_review):
    fh = tapas_physio_get_default_fig_params(verbose);
else
    % Backwards compatibility:
    fh = tapas_physio_get_default_fig_params();
end

if nargin < 4
    physio = tapas_physio_new();
    censoring_threshold = physio.model.movement.censoring_threshold;
end

if nargin < 3
    censoring = [];
end

stringTitle = 'Model: Motion Quality Control - Framewise Displacement';
set(fh, 'Name', stringTitle);

colors = [
    1 0 0
    0 1 0
    0 0 1
    ];

hasOutliers = ~isempty(censoring) && ~isempty(censoring.R_outlier);

% create time vector
nVols = size(rp,1);

t = 1:nVols;
if hasTR
    if numel(TR) > 1 % input is already a time vector
        t = TR;
    else
        t = TR*t - TR;
    end
end

nPlots = 2;
if hasOutliers
    nPlots = 3;
end

%% Realignment parameter
hs(1) = subplot(nPlots,1,1);

for iDim = 1:3
    plot(t, rp(:,iDim), 'Color', colors(iDim,:)); hold all;
    plot(t, quality_measures.rHead*rp(:,iDim+3), 'Color', colors(iDim,:), ...
        'LineStyle', ':');
end
legend('x','pitch','y','roll', 'z', 'yaw');
ylabel('mm');
if ~hasTR
    set(gca,'Xticklabel',[]);
end
title(sprintf('Realignment Parameter (mm), rotation scaled to rHead = %d mm', ...
    quality_measures.rHead));


%% Framewise displacement and friends, subject measures
hs(2) = subplot(nPlots,1,2);
plot(t, quality_measures.absTransDisplacement, 'k'); hold all;
plot(t, quality_measures.absRotDisplacement, 'k--');
plot(t, quality_measures.FD, 'r', 'LineWidth', 4);
plot(t, ones(nVols,1)*censoring_threshold, 'r--')
legend('Absolute Transl. Displacement', 'Absolute Rot. Displacement', ...
    'FD', 'Outlier Threshold')
ylabel('mm');
if ~hasTR
    set(gca,'Xticklabel',[]);
end
title({
    sprintf('Framewise Displacement (mm) and censoring threshold (%.1f)', ...
    censoring_threshold)
    sprintf('Mean FD: %.3f mm; RMS Movement: %.3f', ...
    quality_measures.meanFD, quality_measures.rmsMovement)
    });


%% mask of outlier regressors (stick/spike) for censoring
if hasOutliers
    hs(3) = subplot(nPlots,1,3);
    imagesc(censoring.R_outlier.')
    title('Outlier Mask of Stick (Spike) Regressors for censored volumes');
end

if hasTR
    xlabel('Time (s)')
else
    xlabel('Volume #');
end



tapas_physio_suptitle(stringTitle);

linkaxes(hs, 'x');