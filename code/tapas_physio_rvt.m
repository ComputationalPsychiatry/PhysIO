function [rvt, rpulseMax, rpulseMin] = ...
    tapas_physio_rvt(fr, t)
% computes respiratory volume per time from filtered time series
%
%    [rvt, rpulseMax, rpulseMin] = tapas_physio_rvt(fr, t)
%
%
% The respiratory volume/time is computed for every time point by taking
% the amplitude difference of the closest local maximum and minimum of the
% breathing cycle and dividing by the distance between two subsequent
% breathing maxima
%
% Reference:
%   Birn, R.M., Smith, M.A., Jones, T.B., Bandettini, P.A., 2008. 
%       The respiration response function: The temporal dynamics of 
%       fMRI signal fluctuations related to changes in respiration. 
%       NeuroImage 40, 644?654.
%
% IN
%   fr     filtered respiratory amplitude time series
%   t       vector of time points (seconds) heart rate should be calculated
% OUT
%   rvt         respiratory volume per time vector
%   rpulseMax   vector of maximum inhalation time points
%   rpulseMin   vector of minimum inhalation time points
%
% EXAMPLE
%   [rvt, rpulse] = tapas_physio_rvt(fr, t)
%
%   See also tapas_physio_create_rvt_regressor
%
% Author: Lars Kasper
% Created: 2014-01-20
% Copyright (C) 2013 TNU, Institute for Biomedical Engineering, University of Zurich and ETH Zurich.
%
% This file is part of the physIO toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id: tapas_physio_hr.m 354 2013-12-02 22:21:41Z kasperla $

n = length(t);

dt = t(2)-t(1);
dtBreath = round(2/dt); %in seconds, minimum distance between two breaths

% compute breathing "pulses" (occurence times "rpulse" of max inhalation
% times)
thresh_cardiac = [];
thresh_cardiac.min = .1; 

verbose.level = 2;
verbose.fig_handles = [];

rpulseMax = tapas_physio_get_cardiac_pulses(t, fr, ...
thresh_cardiac,'OXY', dtBreath, verbose);
rpulseMin = tapas_physio_get_cardiac_pulses(t, -fr, ...
thresh_cardiac,'OXY', dtBreath, verbose);

if verbose.level 
    figure('Name', 'Respiratory Volume per Time');
    plot(t,fr, 'g'); hold all;
    stem(rpulseMax, ones(size(rpulseMax)),'b');
    stem(rpulseMin, -ones(size(rpulseMin)),'r');
end

rvt = zeros(size(t));
%for i = 1:n
%
%end
