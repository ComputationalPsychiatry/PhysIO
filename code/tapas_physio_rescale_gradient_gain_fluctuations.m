function G = tapas_physio_rescale_gradient_gain_fluctuations(G, ...
    minStepDistanceSamples, doNormalize)
% Removes infrequent gain changes in gradient time courses (i.e. steps)
%
% G = tapas_physio_rescale_gradient_gain_fluctuations(G, ...
%    minStepDistanceSamples);
%
% IN
%   doNormalize     default: true; if true, every gain interval will have
%                   max 1 after rescaling; 
%                   if false: gain of last interval chosen for all other
% OUT
%
% EXAMPLE
%   tapas_physio_rescale_gradient_gain_fluctuations
%
%   See also
%
% Author: Lars Kasper
% Created: 2015-01-11
% Copyright (C) 2015 TNU, Institute for Biomedical Engineering, University of Zurich and ETH Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id$

% Determine gain fluctuations via steps in sliding-window-maximum
if nargin < 3
    doNormalize = 1;
end
verbose                 = true;
minPeakHeight           = 1000; % TODO: remove this heuristic!
ignoreBoundaryPercent   = 30; % for gain estimation in interval margin 
                             % is ignored in case of a slow change
normFactor              = 1;  % gradients normalized to this value                     
                            
n   = minStepDistanceSamples;
mG  = tapas_physio_maxfilter(abs(G), n);
dmG = diff(mG);

% determine positive and negative steps
[~, idxGainPlus] = findpeaks((dmG), 'minpeakDistance', n, ...
    'minpeakheight', minPeakHeight);
[~, idxGainMinus] = findpeaks(-(dmG), 'minpeakDistance', n, ...
    'minpeakheight', minPeakHeight);

% plus gains refer to max-changes in the future (?)
idxGainPlus     = idxGainPlus + n;

% + 1 because of diff
idxGainMinus    = idxGainMinus + 1;
idxGainSwitch   = [idxGainPlus', idxGainMinus']; 

nGainSwitches = numel(idxGainSwitch);


%% Rescale gradients, if gain switches exist
if nGainSwitches > 0
    
    if verbose
        stringTitle = 'Detected Gradient Gain Fluctuations';
        fh = tapas_physio_get_default_fig_params();
        set(gcf, 'Name', stringTitle);
        hs(1) = subplot(2,1,1);
        hp(1) = plot(G); hold all;
        hp(2) = plot(mG);
        hp(3) = plot(dmG);
        hp(4) = stem(idxGainPlus, mG(idxGainPlus));
        hp(5) = stem(idxGainMinus, mG(idxGainMinus));
    end

    % Sort gain switches and add start/end of gradient time course as
    % interval end points
    idxGainSwitch = [1, sort(idxGainSwitch, 'ascend'), numel(mG)+1];
    
    gainArray = zeros(nGainSwitches+1,1);
    
    % for each gain interval, determine median gain and rescale gradient 
    % time course in interval to gain of last interval
    % last interval, since start of logfile usually not reliably due to
    % dummies, prep-pulses etc.
    for iGainSwitch = (nGainSwitches+1):-1:1
        idxGainStart    = idxGainSwitch(iGainSwitch);
        idxGainEnd      = idxGainSwitch(iGainSwitch+1)-1;
        
        
        ignoreBoundarySamples = ceil((idxGainEnd-idxGainEnd).*...
            ignoreBoundaryPercent/100);
        % Determine gain as max of abs gradient in interval, but ignore
        % transition boundaries
        gainArray(iGainSwitch) = max(abs(G(...
            (idxGainStart+ignoreBoundarySamples): ...
            (idxGainEnd-ignoreBoundarySamples))));
        
         G(idxGainStart:idxGainEnd) = G(idxGainStart:idxGainEnd)/...
             gainArray(iGainSwitch)*normFactor;
         
         if ~doNormalize
             normFactor = gainArray(end);
         end
         
         if verbose
            hp(6) = line([idxGainStart, idxGainEnd], ... 
                 [gainArray(iGainSwitch), gainArray(iGainSwitch)], ...
                 'LineStyle', '--', 'Color',[0 0 0]);
         end
         
    end
    
    if verbose
        legend(hp, {'G','maxFilterG', 'diff maxFilterG', 'Gain Increases', ...
            'Gain Drops', 'Median Gains Per Interval'});
        hs(2) = subplot(2,1,2);

        hp(7) = plot(G, 'LineWidth', 4);
        stringTitle = 'Corrected Gradient Time-course';
        legend(stringTitle);
        title(stringTitle);
        linkaxes(hs, 'x');
    end
  
end