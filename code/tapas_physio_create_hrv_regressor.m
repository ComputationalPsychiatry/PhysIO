function [convHRV, hr, verbose] = tapas_physio_create_hrv_regressor(...
    ons_secs, sqpar, verbose)
% computes cardiac response function regressor and heart rate
%
%    [convHRV, hr] = tapas_physio_create_hrv_regressor(ons_secs, sqpar )
%
% Reference:
%   Chang, Catie, John P. Cunningham, and Gary H. Glover. ???Influence of Heart Rate on the BOLD Signal: The Cardiac Response Function.??? NeuroImage 44, no. 3 (February 1, 2009): 857???869. doi:10.1016/j.neuroimage.2008.09.029.
%
% IN
%   ons_secs.
%       cpulse          onset times (seconds) of heartbeat pulses (R-wave peak)
%       spulse_per_vol  See also tapas_physio_get_sample_points
%   sqpar.
%       onset_slice
% OUT
%   convHRV             cardiac response function regressor after convolution . See
%                       also
% EXAMPLE
%   [convHRV, hr] = tapas_physio_create_hrv_regressor(physio_out.ons_secs, physio_out.sqpar);
%
%   See also tapas_physio_hr tapas_physio_crf
%
% Author: Lars Kasper
% Created: 2013-07-26
% Copyright (C) 2013 TNU, Institute for Biomedical Engineering, University of Zurich and ETH Zurich.
%
% This file is part of the physIO toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id$
if nargin < 3
    verbose.level = [];
    verbose.fig_handles = [];
end

slicenum = 1:sqpar.Nslices;

sample_points  = tapas_physio_get_sample_points(ons_secs, sqpar, slicenum);
hr = tapas_physio_hr(ons_secs.cpulse, sample_points);

if verbose.level>=2
    verbose.fig_handles(end+1) = figure('Name', 'Regressors Heart Rate: HRV X CRF');
    subplot(2,2,1)
    plot(sample_points,hr);xlabel('time (seconds)');ylabel('heart rate (bpm)');
end

% create convolution for whole time series first...
dt = sqpar.TR/sqpar.Nslices;
t = 0:dt:32; % 32 seconds regressor
crf = tapas_physio_crf(t);
crf = crf/max(abs(crf));
% crf = spm_hrf(dt);
if verbose.level>=2
    subplot(2,2,2)
    plot(t, crf);xlabel('time (seconds)');ylabel('cardiac response function');
end

% NOTE: the removal of the mean was implemented to avoid over/undershoots
% at the 1st and last scans of the session due to convolution
convHRV = conv(hr-mean(hr), crf, 'same');

if verbose.level>=2
    subplot(2,2,3)
    plot(sample_points, convHRV);xlabel('time (seconds)');ylabel('heart rate X cardiac response function');
end


% resample to slices needed
hr = hr(sqpar.onset_slice:sqpar.Nslices:end);
convHRV = convHRV(sqpar.onset_slice:sqpar.Nslices:end);
sample_points = sample_points(sqpar.onset_slice:sqpar.Nslices:end);

if verbose.level>=2
    subplot(2,2,4)
    plot(sample_points, convHRV); hold all;
    plot(sample_points, hr);
    xlabel('time (seconds)');ylabel('regessor');
    legend('cardiac response regressor', 'heart rate (bpm)');
end
