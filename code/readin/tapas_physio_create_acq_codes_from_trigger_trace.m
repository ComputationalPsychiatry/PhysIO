function [acq_codes, verbose] = ...
    tapas_physio_create_acq_codes_from_trigger_trace(t, trigger_trace, verbose, ...
    thresholdTrigger, isAlternating)
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
%   verbose         verbosity structure, determines text and graphic output
%                   level, see tapas_physio_new
%   thresholdTrigger    determines what absolute value (difference) should
%                   be consider a trigger onset. default: 1
%                   (jumping from 0 to 1)
%   isAlternating   if true, assume that trigger switches between two levels, e.g.,
%                   +5V and 0 for start of one volume, and back to 5V at start of next volume
%
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

if nargin < 4
    thresholdTrigger = 1;
end

if nargin < 5
    isAlternating = false;
end

isNewTriggerDetection = true;
if isNewTriggerDetection
    % new implementation for ADInstruments, make standard also for BIDS (s.b.
    % in else part); this should also work for non-alternating (see switch
    % below), but does not at the moment for BIDS example

    acq_codes = [];
    iAcqOn = (trigger_trace >= thresholdTrigger); % trigger is 5V, but flips on/off between volumes

    if ~isempty(iAcqOn) % otherwise, nothing to read ...
        % iAcqOn is a column of 1s and 0s, 1 whenever scan acquisition is on
        % Determine 1st start and last stop directly via first/last 1
        % Determine everything else in between via difference (go 1->0 or 0->1)
        iAcqStart   = find(iAcqOn, 1, 'first');
        iAcqEnd     = find(iAcqOn, 1, 'last');
        d_iAcqOn    = diff(iAcqOn);

        % index shift + 1, since diff vector has index of differences i_(n+1) - i_n,
        % and the latter of the two operands (i_(n+1)) has sought value +1
        iAcqStart   = [iAcqStart; find(d_iAcqOn == 1) + 1];
        % no index shift, for the same reason
        iAcqEnd     = [find(d_iAcqOn == -1); iAcqEnd];

        % remove duplicate entries
        iAcqStart = unique(iAcqStart);
        iAcqEnd = unique(iAcqEnd);

        if isAlternating
            % choose all rising and falling triggers as volume starts (instead
            % of interpreting falling as an end of a trigger)
            iAcqStartNew = sort([iAcqStart;iAcqEnd]);
            iAcqEndNew = [];
        else
            iAcqStartNew = iAcqStart;
            iAcqEndNew = iAcqEnd;
        end

        doDebugTriggers = verbose.level >= 3;
        if doDebugTriggers
            verbose.fig_handles(end+1,1) = tapas_physio_get_default_fig_params();
            plot(t,trigger_trace)
            hold all;
            stem(t(iAcqStart), trigger_trace(iAcqStart))
            stem(t(iAcqEnd), trigger_trace(iAcqEnd))
            stem(t(iAcqStartNew), 0.8*trigger_trace(iAcqStartNew))
            strTitle = 'Read-In: Trigger Debugging';
            title(strTitle);
            set(gcf, 'Name', strTitle);
            legend('Volume Trigger Signal', 'Trigger Rising Flank', ...
                'Trigger Falling Flank', 'Chosen Trigger Volume Start')
        end

        iAcqStart = iAcqStartNew;
        iAcqEnd = iAcqEndNew;

        nSamples = size(trigger_trace,1);
        acq_codes = zeros(nSamples,1);
        acq_codes(iAcqStart) = 8; % to match Philips etc. format, volume onset start

        % to match Philips etc. format, volume onset stop,but don't overwrite
        % onsets:
        acq_codes(setdiff(iAcqEnd, iAcqStart)) = 16;

        % report estimated onset gap between last slice of volume_n and 1st slice of
        % volume_(n+1)
        nAcqStarts = numel(iAcqStart);
        nAcqEnds = numel(iAcqEnd);
        nAcqs = min(nAcqStarts, nAcqEnds);

        nDistinctStartEnds = numel(setdiff(iAcqEnd, iAcqStart));
        if nAcqs >= 1 && (nDistinctStartEnds > 0.5*nAcqStarts)
            % report time of acquisition, as defined in SPM, but only if
            % majority of the acquisition trigger starts don't coincide with the
            % end triggers (i.e., length 1 triggers)
            TA = mean(t(iAcqEnd(1:nAcqs)) - t(iAcqStart(1:nAcqs)));
            verbose = tapas_physio_log(...
                sprintf('TA = %.4f s (Estimated time of acquisition during one volume TR)', ...
                TA), verbose, 0);
        end
    end

else
    % OLD implementation from BIDS, for some reason the acq_codes 16 above change the reference results
    % TODO: unify both trigger computations

    acq_codes = [];
    iAcqStart = (trigger_trace >= thresholdTrigger); % trigger has 1, rest is 0;
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


         doDebugTriggers = verbose.level >= 3;
        if doDebugTriggers
            verbose.fig_handles(end+1,1) = tapas_physio_get_default_fig_params();
            plot(t,trigger_trace)
            hold all;
            stem(t(iAcqStart), trigger_trace(iAcqStart))
            strTitle = 'Read-In: Trigger Debugging';
            title(strTitle);
            set(gcf, 'Name', strTitle);
            legend('Volume Trigger Signal', 'Trigger Rising Flank')
        end

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
end