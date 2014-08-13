function [cpulse, verbose] = tapas_physio_findpeaks_template_correlation(...
    c, pulseCleanedTemplate, cpulseSecondGuess, averageHeartRateInSamples, ...
    verbose)
% Finds peaks of a time series via pre-determined template via maxima of
% correlations via going backward from search starting point in time
% series, and afterwards forward again
%
%   [cpulse, verbose] = tapas_physio_findpeaks_template_correlation(...
%       c, pulseCleanedTemplate, cpulseSecondGuess, averageHeartRateInSamples, verbose)
%
% IN
%
% OUT
%
% EXAMPLE
%   tapas_physio_findpeaks_template_correlation
%
%   See also
%
% Author: Steffen Bollmann, merged in this function: Lars Kasper
% Created: 2014-08-05
% Copyright (C) 2014 TNU, Institute for Biomedical Engineering, University of Zurich and ETH Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id$

% Determine starting peak for the search:
%   search for a representative R-peak in the first 20 peaks

% TODO: maybe replace via convolution with template? "matched
% filter theorem"?
debug = verbose.level >= 4;

halfTemplateWidthInSamples = floor(numel(pulseCleanedTemplate)/2);

centreSampleStart = round(2*halfTemplateWidthInSamples+1);
centreSampleEnd = cpulseSecondGuess(20);
for n=centreSampleStart:centreSampleEnd
    startSignalIndex=n-halfTemplateWidthInSamples;
    endSignalIndex=n+halfTemplateWidthInSamples;
    
    signalPart = c(startSignalIndex:endSignalIndex);
    correlation = corrcoef(signalPart,pulseCleanedTemplate);
    
    %Debug
    if debug && ~mod(n, 100)
        figure(2);clf;
        plot(signalPart);
        hold all;
        plot(pulseCleanedTemplate);
    end
    % Debug
    
    similarityToTemplate(n) = correlation(1,2);
end

[C_bestMatch,I_bestMatch] = max(similarityToTemplate);
clear similarityToTemplate

%% now compute backwards to the beginning:
% go average heartbeat by heartbeat back and look (with
% decreasing weighting for higher distance) for highest
% correlation with template heartbeat

n=I_bestMatch;
bestPosition = n; % to capture case where 1st R-peak is best

peakNumber = 1;
similarityToTemplate=zeros(size(c,1),1);

searchStepsTotal=round(0.5*averageHeartRateInSamples);
while n > 1+searchStepsTotal+halfTemplateWidthInSamples
    for searchPosition=-searchStepsTotal:1:searchStepsTotal
        startSignalIndex=n-halfTemplateWidthInSamples+searchPosition;
        endSignalIndex=n+halfTemplateWidthInSamples+searchPosition;
        
        signalPart = c(startSignalIndex:endSignalIndex);
        correlation = corrcoef(signalPart,pulseCleanedTemplate);
        
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
        correlationWeighted =  currentWeight .* correlation(1,2);
        similarityToTemplate(n+searchPosition) = correlationWeighted;
        
    end
    
    %DEBUG
    if debug
        if (n>100) && (n< size(c,1)-100)
            figure(1);subplot 211;plot(t,similarityToTemplate,'b-')
            title('Finding first peak (heartbeat/max inhale), backwards');
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
    
    cpulse(peakNumber) = bestPosition;
    peakNumber = peakNumber+1;
    
    
    n=bestPosition-averageHeartRateInSamples;
end % END: going backwards to beginning of time course

%% Now go forward through the whole time series
n=bestPosition; % 1st R-peak
peakNumber=1;
clear cpulse;
%now correlate template with PPU signal at the positions
%where we would expect a peak based on the average heartrate and
%search in the neighborhood for the best peak, but weight the peaks
%deviating from the initial starting point by a gaussian
searchStepsTotal=round(0.5*averageHeartRateInSamples);

if n < searchStepsTotal+halfTemplateWidthInSamples+1
    n=searchStepsTotal+halfTemplateWidthInSamples+1;
end

while n < size(c,1)-searchStepsTotal-halfTemplateWidthInSamples
    %search around peak
    
    for searchPosition=-searchStepsTotal:1:searchStepsTotal
        startSignalIndex=n-halfTemplateWidthInSamples+searchPosition;
        endSignalIndex=n+halfTemplateWidthInSamples+searchPosition;
        
        signalPart = c(startSignalIndex:endSignalIndex);
        correlation = corrcoef(signalPart,pulseCleanedTemplate);
        
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
        
        gaussianWindow = gausswin(2*searchStepsTotal+1);
        locationWeight = gaussianWindow(searchPosition+searchStepsTotal+1);
        %                     locationWeight = 1;
        amplitudeWeight = abs(c(n+searchPosition+1));
        %                     amplitudeWeight = 1;
        correlationWeighted =  locationWeight .* amplitudeWeight .* correlation(1,2);
        similarityToTemplate(n+searchPosition) = correlationWeighted;
        
        
    end
    
    %DEBUG
    if debug
        if (n>100) && (n< size(c,1)-100)
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
    
    if  (foundCpulses < 21) && (foundCpulses >= 3)
        currentHeartRateInSamples = round(mean(diff(cpulse)));
    end
    
    if foundCpulses > 20
        currentCpulses = cpulse (foundCpulses-20:foundCpulses);
        currentHeartRateInSamples = round(mean(diff(currentCpulses)));
    end
    
    
    %check currentHeartRate
    checkSmaller=currentHeartRateInSamples > 0.5*averageHeartRateInSamples;
    checkLarger=currentHeartRateInSamples < 1.5*averageHeartRateInSamples;
    
    %jumpToNextPeakSearchArea
    if (checkSmaller && checkLarger)
        n=bestPosition+currentHeartRateInSamples;
    else
        n=bestPosition+averageHeartRateInSamples;
    end
end