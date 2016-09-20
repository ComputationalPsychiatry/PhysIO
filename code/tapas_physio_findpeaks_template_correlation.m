function [cpulse, verbose, plotData] = tapas_physio_findpeaks_template_correlation(...
    c, pulseCleanedTemplate, cpulseSecondGuess, averageHeartRateInSamples, ...
    verbose, varargin)
% Finds peaks of a time series via pre-determined template via maxima of
% correlations via going backward from search starting point in time
% series, and afterwards forward again
%
%   [cpulse, verbose, plotData] = tapas_physio_findpeaks_template_correlation(...
%       c, pulseCleanedTemplate, cpulseSecondGuess, averageHeartRateInSamples, verbose)
%
% IN
%   varargin    property name/value pairs for additional options
%
%   'idxStartPeakSearch' [1,2] indices of cpulseSecondGuess 
%                       giving start and end point of search range for
%                       detection of representative starting cycle, from
%                       which both cycles before and after will be
%                       determined via the pulseCleanedTemplate
%
% OUT
%   cpulse
%   verbose
%   plotData        structure of plot-relevant data for algorithm depiction
%                   fields: searchedAt, amplitudeWeight, locationWeight,
%                   similarityToTemplate
%
% EXAMPLE
%   tapas_physio_findpeaks_template_correlation
%
%   See also
%
% Author: Steffen Bollmann, merged in this function: Lars Kasper
% Created: 2014-08-05
% Copyright (C) 2014 TNU, Institute for Biomedical Engineering, 
%                         University of Zurich and ETH Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id$

% Determine starting peak for the search:
%   search for a representative R-peak a range of peaks

% TODO: maybe replace via convolution with template? "matched
% filter theorem"?


nSamples = size(c,1);

debug = verbose.level >= 4;

defaults.idxStartPeakSearch = [0 20];
defaults.t = 1:nSamples;

args = tapas_physio_propval(varargin, defaults);
tapas_physio_strip_fields(args);

halfTemplateWidthInSamples = floor(numel(pulseCleanedTemplate)/2);

[tmp,zTransformedTemplate] = tapas_physio_corrcoef12(pulseCleanedTemplate,...
    pulseCleanedTemplate);
isZTransformed = [0 1];

%% Find best (representative) R-peak within first 20(=idxStartPeakSearch) 
% cycles to start backwards search
% start and end point of search for representative start cycle
centreSampleStart = round(2*halfTemplateWidthInSamples+1);

if idxStartPeakSearch(1) > 0
    centreSampleStart = centreSampleStart + ...
        cpulseSecondGuess(idxStartPeakSearch(1));
end
centreSampleEnd = cpulseSecondGuess(idxStartPeakSearch(2));


similarityToTemplate = zeros(1, ceil(centreSampleEnd));
for n=centreSampleStart:centreSampleEnd
    startSignalIndex    = n - halfTemplateWidthInSamples;
    endSignalIndex      = n + halfTemplateWidthInSamples;
    
    signalPart = c(startSignalIndex:endSignalIndex);
    similarityToTemplate(n) = tapas_physio_corrcoef12(signalPart, ...
        zTransformedTemplate, isZTransformed);
    
    %Debug
    if debug && ~mod(n, 100)
        figure(2);clf;
        plot(signalPart);
        hold all;
        plot(pulseCleanedTemplate);
    end
    % Debug
    
end

[C_bestMatch, I_bestMatch] = max(similarityToTemplate);
clear similarityToTemplate



%% now compute backwards to the beginning:
% go average heartbeat by heartbeat back and look (with
% decreasing weighting for higher distance) for highest
% correlation with template heartbeat

plotData.searchedAt = zeros(size(c));
plotData.locationWeight  = zeros(size(c));
plotData.amplitudeWeight  = zeros(size(c));

n = I_bestMatch;
bestPosition = n; % to capture case where 1st R-peak is best

peakNumber = 1;
similarityToTemplate = zeros(nSamples,1);

searchStepsTotal = round(0.5*averageHeartRateInSamples);

% zero-pad phys time course to allow computing correlations for first heartbeat 
cPadded = [zeros(halfTemplateWidthInSamples,1); c];
bestPosition = bestPosition + halfTemplateWidthInSamples; % index update for padded c

% Stepping backwards through heartbeat intervals
while n > 1+searchStepsTotal+halfTemplateWidthInSamples
    
    % search samples within current heartbeat interval
    for searchPosition = -searchStepsTotal:1:searchStepsTotal
        startSignalIndex    = n - halfTemplateWidthInSamples+searchPosition;
        endSignalIndex      = n + halfTemplateWidthInSamples+searchPosition;
        
        signalPart          = cPadded(startSignalIndex:endSignalIndex);
        correlation = tapas_physio_corrcoef12(signalPart,zTransformedTemplate, ...
            isZTransformed);
        
        %DEBUG
        if debug && ~mod(n, 100)
            figure(1);
            subplot 212;
            plot(signalPart);
            hold all;
            plot(pulseCleanedTemplate);
            hold off;
            title('Correlating current window with template wave');
        end
        %DEBUG
        
        % weight correlations far away from template center
        % less; since heartbeat to be expected in window center
        % gaussianWindow = gausswin(2*searchStepsTotal+1);
        %                     currentWeight = gaussianWindow(searchPosition+searchStepsTotal+1);
        
        currentWeight = abs(c(n+searchPosition+1));
        correlationWeighted =  currentWeight .* correlation;
        similarityToTemplate(n+searchPosition) = correlationWeighted;
        
    end % search within heartbeat
    
    %DEBUG
    if debug
        if (n>100) && (n< nSamples-100)
            figure(1);subplot 211;plot(t,similarityToTemplate,'b-');
            title('Finding first peak (heartbeat/max inhale), backwards');
            xlim([t(n-100) t(n+100)]);
        end
    end
    %DEBUG
    
    
    %find biggest correlation-peak from the last search
    indexSearchStart=n-searchStepsTotal;
    indexSearchEnd=n+searchStepsTotal;
    
    indexSearchRange=indexSearchStart:indexSearchEnd;
    searchRangeValues=similarityToTemplate(indexSearchRange);
    [C_bestMatch,I_bestMatch] = max(searchRangeValues);
    bestPosition = indexSearchRange(I_bestMatch); 
    
    n=bestPosition-averageHeartRateInSamples;
end % END: going backwards to beginning of time course

%% Now go forward through the whole time series
% 1st R-peak, correct for index change with zero-padding
n           = bestPosition - halfTemplateWidthInSamples; 
peakNumber  = 1;
clear cpulse;

% Now correlate template with PPU signal at the positions
% where we would expect a peak based on the average heartrate and
% search in the neighborhood for the best peak, but weight the peaks
% deviating from the initial starting point by a gaussian
searchStepsTotal = round(0.5*averageHeartRateInSamples);

% for weighted searching of max correlation
gaussianWindow = gausswin(2*searchStepsTotal+1);

n = max(n, searchStepsTotal + halfTemplateWidthInSamples + 1);

% zero-pad c at end to allow for detection of last peak by
% template-matching up to the last sample of c
c = [c; zeros(searchStepsTotal + halfTemplateWidthInSamples + 1, 1)];

while n < nSamples % -searchStepsTotal - halfTemplateWidthInSamples
    
    %search around peak
    for searchPosition = -searchStepsTotal:1:searchStepsTotal
        startSignalIndex = ...
            max(1, n-halfTemplateWidthInSamples+searchPosition);
        endSignalIndex = n+halfTemplateWidthInSamples+searchPosition;
        
        signalPart = c(startSignalIndex:endSignalIndex);
            
        
        correlation = tapas_physio_corrcoef12(signalPart, ...
            zTransformedTemplate, isZTransformed);
        
        %DEBUG
        if debug && ~mod(n, 100)
            figure(1);
            subplot 212;
            plot(signalPart);
            hold all;
            plot(pulseCleanedTemplate);
            hold off;
        end
        %  DEBUG
        
        locationWeight = gaussianWindow(searchPosition+searchStepsTotal+1);
        %                     locationWeight = 1;
        amplitudeWeight = abs(c(n+searchPosition+1));
        %                     amplitudeWeight = 1;
        correlationWeighted =  locationWeight .* amplitudeWeight .* correlation;
        similarityToTemplate(n+searchPosition) = correlationWeighted;
        
        % collect plot Data
        plotData.locationWeight(n+searchPosition) = locationWeight;
        plotData.amplitudeWeight(n+searchPosition) = amplitudeWeight;
        
        if searchPosition==0
            plotData.searchedAt(n+searchPosition) = 1;
        end
        
    end
    
    %DEBUG
    if debug
        if (n>100) && (n< nSamples-100)
            figure(1);subplot 211;plot(t,similarityToTemplate,'b-')
            xlim([t(n-100) t(n+100)])
        end
    end
    %DEBUG
    
    
    %find biggest correlation-peak from the last search
    indexSearchStart=n-searchStepsTotal;
    indexSearchEnd=n+searchStepsTotal;
    
    indexSearchRange=indexSearchStart:indexSearchEnd;
    searchRangeValues=similarityToTemplate(indexSearchRange);
    [C_bestMatch,I_bestMatch] = max(searchRangeValues);
    bestPosition = indexSearchRange(I_bestMatch);
    
    if debug
        stem(t(bestPosition),4,'g');
    end
    
    cpulse(peakNumber) = bestPosition;
    peakNumber = peakNumber+1;
    
    %only take the last 20 cpulses to compute the current HeartRate
    foundCpulses = size(cpulse,2);
    
    if  foundCpulses < 3
        currentHeartRateInSamples=averageHeartRateInSamples;
    end
    
    if (foundCpulses < 21) && (foundCpulses >= 3)
        currentHeartRateInSamples = round(mean(diff(cpulse)));
    end
    
    if foundCpulses >= 21
        currentCpulses = cpulse (foundCpulses-20:foundCpulses);
        currentHeartRateInSamples = round(mean(diff(currentCpulses)));
    end
    
    
    %check currentHeartRate
    checkSmaller    = currentHeartRateInSamples > 0.5*averageHeartRateInSamples;
    checkLarger     = currentHeartRateInSamples < 1.5*averageHeartRateInSamples;
    
    %jumpToNextPeakSearchArea
    if (checkSmaller && checkLarger)
        n=bestPosition+currentHeartRateInSamples;
    else
        n=bestPosition+averageHeartRateInSamples;
    end
end

plotData.similarityToTemplate = similarityToTemplate;
