function [rvt] = tapas_physio_rvt_hilbert(fr, t, sample_points, verbose)
% computes respiratory volume per time from filtered time series
%
%    [rvt] = tapas_physio_rvt(fr, t)
%
% The respiratory volume/time is computed by calculating the instantaneous
% amplitude / frequency of the breathing signal via the Hilbert transform.
%
% Reference:
%   Birn, R.M., Smith, M.A., Jones, T.B., Bandettini, P.A., 2008.
%       The respiration response function: The temporal dynamics of
%       fMRI signal fluctuations related to changes in respiration.
%       NeuroImage 40, 644-654.
%
% IN
%   fr     filtered respiratory amplitude time series
%   t      time vector for fr
%   sample_points       vector of time points (seconds) respiratory volume/time should be calculated
% OUT
%   rvt         respiratory volume per unit time vector
%
% EXAMPLE
%   [rvt, rpulse] = tapas_physio_rvt(fr, t)
%
%   See also tapas_physio_create_rvt_regressor

% Author: Sam Harrison
% Created: 2019-05-10
% Copyright (C) 2019 TNU, Institute for Biomedical Engineering, University of Zurich and ETH Zurich.
%
% This file is part of the physIO toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

if nargin < 3
    sample_points = t;
end
if nargin < 4
    verbose.level = 0;
    verbose.fig_handles = [];
end

f_sample = 1 / (t(2)-t(1));
n_pad = ceil(10.0 * f_sample);

%% Respiratory volume is amplitude envelope %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mild low-pass filter to remove high-frequency noise
d = designfilt( ...
    'lowpassiir', 'FilterOrder', 10, ...
    'HalfPowerFrequency', 2.0, 'SampleRate', f_sample);
fr_lp = filtfilt(d, padarray(fr, n_pad, 'circular'));
fr_lp = fr_lp(n_pad+1:end-n_pad);

% Analytic signal -> magnitude
fr_analytic = hilbert(fr_lp);
fr_mag = abs(fr_analytic);

% Low-pass filter envelope to retrieve change in respiratory volume
d = designfilt( ...
    'lowpassiir', 'FilterOrder', 10, ...
    'HalfPowerFrequency', 0.2, 'SampleRate', f_sample);
fr_rv = filtfilt(d, padarray(fr_mag, n_pad, 'circular'));
fr_rv = fr_rv(n_pad+1:end-n_pad);
fr_rv(fr_rv < 0.0) = 0.0;

if verbose.level>=2
    verbose.fig_handles(end+1) = tapas_physio_get_default_fig_params();
    set(gcf, 'Name', 'Model: Respiratory Volume');
    hold all;
    hp(1) = plot(t, fr);
    hp(2) = plot(t, fr_lp);
    hp(3) = plot(t, fr_mag);
    hp(4) = plot(t, fr_rv);
    strLegend = {
        'Filtered breathing signal', ...
        '... after low pass-filter', ...
        'Breathing signal envelope', ...
        'Respiratory volume'};
    legend(hp, strLegend)
end

%% Breathing rate is instantaneous frequency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Low-pass filter at approximately breathing-rate
d = designfilt( ...
    'lowpassiir', 'FilterOrder', 10, ...
    'HalfPowerFrequency', 0.75, 'SampleRate', f_sample);
fr_lp = filtfilt(d, padarray(fr, n_pad, 'circular'));
fr_lp = fr_lp(n_pad+1:end-n_pad);

% Now iteratively refine instantaneous frequency estimate
% Aim is to remove any high frequencies caused by funny shaped waveforms
fr_filt = fr_lp;
for n = 1:5
    % Analytic signal -> phase
    fr_analytic = hilbert(fr_filt);
    fr_phase = phase(fr_analytic);
    
    % Remove any phase decreases that may occur
    % Find places where the gradient changes sign
    fr_phase_diff = diff(sign(gradient(fr_phase)));
    decrease_inds = find(fr_phase_diff < 0);
    increase_inds = [find(fr_phase_diff > 0); length(fr_phase)];
    for n_start = decrease_inds'
        %   /2\   /4
        % 1/   \3/
        % Find value of `fr` at:
        %   [2]: start (i.e. peak)
        %   [3]: end (i.e. trough)
        fr_start = fr_phase(n_start);
        n_end = increase_inds(find(increase_inds > n_start, 1));
        fr_end = fr_phase(n_end);
        
        % Now find where `fr` passes `fr_end` for the first time [1]
        n_min = find(fr_phase > fr_end, 1);
        if isempty(n_min)
            n_min = n_start;
        end
        % And find where `fr` passes `fr_end` for the last time [4]
        n_max = find(fr_phase < fr_start, 1, 'last');
        if isempty(n_max)
            n_max = n_end;
        end
        
        % Finally, linearly interpolate from [1] to [4]
        fr_phase(n_min:n_max) = linspace(fr_end, fr_start, n_max-n_min+1);
    end
    
    % And filter out any high frequencies from phase-only signal
    fr_filt = filtfilt(d, padarray(cos(fr_phase), n_pad, 'circular'));
    fr_filt = fr_filt(n_pad+1:end-n_pad);
end

% Recalculate analytic signal -> phase
fr_phase = phase(hilbert(fr_filt));

% Transform to instantaneous frequency
fr_if = f_sample * gradient(fr_phase) / (2 * pi);
fr_if(fr_if >  2.0) = 2.0;   % Lower limit of 2.0 breaths per second
fr_if(fr_if < 0.05) = 0.05;  % Upper limit of 20.0 s per breath

if verbose.level>=2
    verbose.fig_handles(end+1) = tapas_physio_get_default_fig_params();
    set(gcf, 'Name', 'Model: Breathing rate');
    hold all;
    hp(1) = plot(t, fr);
    hp(2) = plot(t, fr_lp);
    hp(3) = plot(t, std(fr) * cos(fr_phase));
    hp(4) = plot(t, fr_if);
    strLegend = {
        'Filtered breathing signal', ...
        '... after low pass-filter', ...
        '... after removing amplitude', ...
        'Instantaneous breathing rate'};
    legend(hp, strLegend)
end

% figure; hold all; plot(t, fr); plot(t, fr_rv .* cos(fr_phase)); plot(t, fr_mag .* cos(fr_phase));
% figure; hold all; plot(abs(fft(zscore(fr_rv)))); plot(abs(fft(zscore(fr_if))));

%% And make RVT! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RVT = magnitude * breathing rate
fr_rvt = fr_rv .* fr_if;

% Need to downsample to `sample_points`, taking care to avoid aliasing
f_sample_out = 1 / mean(diff(sample_points));
[re_rvt, t_re_rvt] = resample(fr_rvt, t, f_sample_out);
% And now interpolate onto new timepoints
rvt = interp1(t_re_rvt, re_rvt, sample_points, 'linear');
% Nearest-neighbour interpolation for outside recorded timepoints
% Be more careful here as don't want RVT to go negative
if sum(isnan(rvt)) > 0
    nan_inds = isnan(rvt);
    rvt(nan_inds) = interp1( ...
        t_re_rvt, re_rvt, sample_points(nan_inds), ...
        'nearest', 'extrap');
end

% figure; hold all; plot(t, fr_rvt); plot(sample_points, rvt);

end