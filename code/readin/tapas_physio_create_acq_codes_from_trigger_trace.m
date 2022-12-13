function [acq_codes, verbose] = ...
    tapas_physio_create_acq_codes_from_trigger_trace(t, trigger_trace, verbose)
% Creates integer acquisition codes (on/off volume/slice start/end events)
% from continuous trigger trace (e.g., TTL trigger spiking from 0 to 5V,
% or alternating between two different voltage levels)
%
% [acq_codes, verbose] = ...
%     tapas_physio_create_acq_codes_from_trigger_trace(t, trigger_trace, verbose)
%
% IN
%   t               [nSamples,1] time vector corresponding to trigger trace
%   trigger_trace   [nSamples, 1] trigger trace (e.g., TTL 0 to 5 V)
% OUT
%
% EXAMPLE
%   tapas_physio_create_acq_codes_from_trigger_trace
%
%   See also

% Author:   Lars Kasper
% Created:  2022-12-13
% Copyright (C) 2022 TNU, Institute for Biomedical Engineering,
%                    University of Zurich and ETH Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under
% the terms of the GNU General Public License (GPL), version 3. You can
% redistribute it and/or modify it under the terms of the GPL (either
% version 3 or, at your option, any later version). For further details,
% see the file COPYING or <http://www.gnu.org/licenses/>.

acq_codes = [];
iAcqStart = (trigger_trace~=0); % trigger has 1, rest is 0;
if ~isempty(iAcqStart) % otherwise, nothing to read ...
    % iAcqStart is a columns of 0 and 1, 1 for the trigger event of a new
    % volume start

    % sometimes, trigger is on for several samples; ignore these extended
    % phases of "on-triggers" as duplicate values, if trigger distance is
    % very different
    %
    % fraction of mean trigger distance; if trigger time difference below that, will be removed
    outlierThreshold = 0.2;

    idxAcqStart = find(iAcqStart);
    dAcqStart = diff(idxAcqStart);

    % + 1 because of diff
    iAcqOutlier = 1 + find(dAcqStart < outlierThreshold*mean(dAcqStart));
    iAcqStart(idxAcqStart(iAcqOutlier)) = 0;

    nSamples = size(trigger_trace,1);
    acq_codes = zeros(nSamples,1);
    acq_codes(iAcqStart) = 8; % to match Philips etc. format

    nAcqs = numel(find(iAcqStart));

    if nAcqs >= 1
        % report time of acquisition, as defined in SPM
        meanTR = mean(diff(t(iAcqStart)));
        stdTR = std(diff(t(iAcqStart)));
        verbose = tapas_physio_log(...
            sprintf('TR = %.3f +/- %.3f s (Estimated mean +/- std time of repetition for one volume; nTriggers = %d)', ...
            meanTR, stdTR, nAcqs), verbose, 0);
    end
end

