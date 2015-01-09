function [VOLLOCS, LOCS, verbose] = ...
    tapas_physio_create_scan_timing_from_gradients_auto_philips(log_files, ...
    sqpar, verbose)
% Automatically extracts slice/volume scan events from gradient timecourse
% SCANPHYSLOG file
%
% [VOLLOCS, LOCS, verbose] = ...
%     tapas_physio_create_scan_timing_from_gradients_auto_philips(log_files, ...
%     sqpar, verbose)
%
% This function determines slice/volume starts from the gradient time course 
% automatically by assuming a regularity of them from the 
% sequence timing parameters, in particular TR and number of slices
%
% Therefore, unlike tapas_physio_create_scan_timing_from_gradients_philips, 
% no thresholds of slice/volume starts have to be given to determine the
% timing and are inferred on iteratively instead.
%
%   Workflow
%   1. Determine template of each volume gradient time course by using TR. 
%   2. Determine volume events using volume template (counting either from
%   start or end of the time series)
%   3. Determine slice events between all detected volumes
%       (again, creating a slice template and matching it to time series
%       between consecutive volumes...)
%
% IN
%   log_files   is a structure containing the following filenames (with full
%           path)
%       .log_cardiac        contains ECG or pulse oximeter time course
%                           for Philips: 'SCANPHYSLOG<DATE&TIME>.log';
%                           can be found on scanner in G:/log/scanphyslog-
%                           directory, one file is created per scan, make sure to take
%                           the one with the time stamp corresponding to your PAR/REC
%                           files
%       .log_respiration    contains breathing belt amplitude time course
%                           for Philips: same as .log_cardiac
%
%
%   sqpar                   - sequence timing parameters
%           .nSlices        - number of slices per volume in fMRI scan
%           .nSlicesPerBeat - usually equals nSlices, unless you trigger with the heart beat
%           .TR             - repetition time in seconds
%           .nDummies       - number of dummy volumes
%           .nScans         - number of full volumes saved (volumes in nifti file,
%                             usually rows in your design matrix)
%           .Nprep          - number of non-dummy, volume like preparation pulses
%                             before 1st dummy scan. If set, logfile is read from beginning,
%                             otherwise volumes are counted from last detected volume in the logfile
%           .time_slice_to_slice - time between the acquisition of 2 subsequent
%                             slices; typically TR/nSlices or
%                             minTR/nSlices, if minimal temporal slice
%                             spacing was chosen
%            onset_slice    - slice whose scan onset determines the adjustment of the
%                             regressor timing to a particular slice for the whole volume
%
%                             NOTE: only necessary, if thresh.grad_direction is empty
%   verbose
%
% OUT
%           VOLLOCS         - locations in time vector, when volume scan
%                             events started
%           LOCS            - locations in time vector, when slice or volume scan
%                             events started
%
% EXAMPLE
%   [VOLLOCS, LOCS] = tapas_physio_create_scan_timing_from_gradients_philips(logfile,
%   thresh.scan_timing);
%
%   See also tapas_physio_create_scan_timing_from_gradients_philips
%
% Author: Lars Kasper
% Created: 2015-01-09
% Copyright (C) 2013 Institute for Biomedical Engineering, ETH/Uni Zurich.
%
% This file is part of the PhysIO toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id: tapas_physio_create_scan_timing_from_gradients_philips.m 632 2015-01-09 12:36:12Z kasperla $


% smaller than typical single shot EPI slice duration (including waiting
% for TE)
minSliceDuration = 0.040; 

doCountSliceEventsFromLogfileStart  = isfield(sqpar, 'Nprep') && ...
    ~isempty(sqpar.Nprep);


% everything stored in 1 logfile
if ~isfield(log_files, 'cardiac') || isempty(log_files.cardiac)
    logfile = log_files.respiration;
else
    logfile = log_files.cardiac;
end

nScans          = sqpar.Nscans;
nDummies        = sqpar.Ndummies;
nSlicesPerBeat  = sqpar.NslicesPerBeat;
nSlices         = sqpar.Nslices;


if doCountSliceEventsFromLogfileStart
    Nprep = sqpar.Nprep;
end

%% Read columns of physlog-file and convert into double
% use textread as long as it exists, for it is much faster (factor 4) than
% textscan; TODO: use fread ans sscanf to make it even faster...
if exist('textread')
    [z{1:10}]   = textread(logfile,'%d %d %d %d %d %d %d %d %d %s','commentstyle', 'shell');
else
    fid     = fopen(logfile, 'r');
    z       = textscan(fid, '%d %d %d %d %d %d %d %d %d %s', 'commentstyle', '#');
    z(1:9)  = cellfun(@double, z(1:9), 'UniformOutput', false);
    fclose(fid);
end

z{10}       = hex2dec(z{10}); % hexadecimal acquisition codes converted;
y           = cell2mat(z);

acq_codes   = y(:,10);
nSamples    = size(y,1);

dt          = log_files.sampling_interval(1);

%default: 500 Hz sampling frequency
if isempty(dt)
    dt      = 2e-3;
end

t           = -log_files.relative_start_acquisition + ((0:(nSamples-1))*dt)';



% finding scan volume starts (svolpulse)
thresh.grad_direction = 'abs';
switch lower(thresh.grad_direction)
    case 'x'
        gradient_choice = y(:,7);
    case 'y'
        gradient_choice = y(:,8);
    case 'z'
        gradient_choice = y(:,9);
    case {'xyz', 'abs'}
        gradient_choice = sqrt(sum(y(:,7:9).^2,2));
end
gradient_choice         = reshape(gradient_choice, [] ,1);


[tmp, LOCS]    = tapas_physio_findpeaks(z2,'minpeakheight',thresh.slice, ...
    'minpeakdistance', ceil(minSliceDuration/dt));



%% 1. Determine template for a gradient time-course during a volume

thresh_min                  = tapas_physio_prctile(gradient_choice, 20);
minVolumeDistanceSamples    = ceil(sqpar.TR*0.95/dt);
[templateGradientVolume, secondGuessVOLLOCS, averageTR] = ...
    tapas_physio_get_cardiac_pulse_template(t, gradient_choice, thresh_min, ...
    minVolumeDistanceSamples, verbose);



%% 2. Determine volume events from template using cross-correlation
 
[VOLLOCS, verbose] = tapas_physio_findpeaks_template_correlation(...
            c, templateGradientVolume, secondGuessVOLLOCS, averageHeartRateInSamples, ...
            verbose);

        
        
%% 3. Determine slice events from volume positions and info on number of slices
        
nVolumes = numel(VOLLOCS);

% Start searching for slice event slightly before volume event to include
% the volume event as the first slice
nShiftSamples = ceil(minSliceDuration/2/dt); 
fprintf('Finding slice events of volumes %04d/%04d',0, nVolumes);
LOCS = cell(nVolumes,1);

minSliceDistanceSamples = ceil(minSliceDuration/dt);
for iVol = 1:nVolumes
    fprintf('\b\b\b\b\b\b\b\b\b%04d/%04d', iVol, nVolumes);
    
    iStart          = VOLLOCS(iVol) - nShiftSamples;
    iEnd            = VOLLOCS(iVol+1) - nShiftSamples;
    [~, LOCS{iVol}] = tapas_physio_findpeaks(...
        gradient_choice(iStart:iEnd), ...
        'minPeakDistance', minSliceDistanceSamples, ...
        'nPeaks', nSlices);

    LOCS{iVol} = LOCS{iVol} + iStart - 1;
end
fprintf('\n');

LOCS = cell2mat(LOCS')';



%% Select relevant events from detected ones using sequence parameter info

try
    
    if doCountSliceEventsFromLogfileStart
        VOLLOCS = LOCS(Nprep*nSlices + ...
            (1:nSlices:(nDummies+nScans)*nSlices));
    else % count from end
        VOLLOCS = LOCS((end-(nDummies+nScans)*nSlices+1):nSlices:end);
    end

    LOCS    = reshape(LOCS, [], 1);
    VOLLOCS = reshape(VOLLOCS, [], 1);
catch
    VOLLOCS = [];
end

if verbose.level>=1
    
    % Depict all gradients, raw
    verbose.fig_handles(end+1) = tapas_physio_get_default_fig_params();
    set(gcf,'Name', 'Thresholding Gradient for slice acq start detection');
    fs(1) = subplot(3,1,1); 
    
    plot(t, sqrt(sum(y(:,7:9).^2,2)), '--k');
    hold all;
    plot(t, y(:,7:9));
    
    
    if ismember(8,acq_codes)
        hold all;
        stem(t, acq_codes*max(max(abs(y(:,7:9))))/20);
    end
    
    
    legend('abs(G_x^2+G_y^2+G_z^2)', 'gradient x', 'gradient y', 'gradient z');
    title('Raw Gradient Time-courses');
    
    % Plot gradient thresholding for slice timing determination
    fs(2) = subplot(3,1,2);
    hp = plot(t,[gradient_choice z2]); hold all;
    hp(end+1) = plot(t, repmat(thresh.zero, nSamples, 1));
    hp(end+1) = plot(t, repmat(thresh.slice, nSamples, 1));
    lg = {'Chosen gradient for thresholding', ...
        'Gradient with values < thresh.zero set to 0', ...
        'thresh.zero', 'thresh.slice'};
    
    if doDetectVolumesByGradientAmplitude
        hp(end+1) = plot(t, repmat(thresh.vol, nSamples, 1));
        lg{end+1} = 'thresh.vol';
    end
    title({'Thresholding Gradient for slice acq start detection', '- found scan events -'});
    legend(hp, lg);
    xlabel('t(s)');
    
    % Plot gradient thresholding for slice timing determination
    
    if ~isempty(VOLLOCS)
        hp(end+1) = stem(t(VOLLOCS), 1.25*max(gradient_choice)*ones(size(VOLLOCS))); hold all
        lg{end+1} = sprintf('Found volume events (N = %d)', numel(VOLLOCS));
    end
    
    if ~isempty(LOCS)
        hp(end+1) = stem(t(LOCS), max(gradient_choice)*ones(size(LOCS))); hold all
        lg{end+1} = sprintf('Found slice events (N = %d)', numel(LOCS));
        
        dLocsSecs = diff(LOCS)*dt*1000;
        ymin = tapas_physio_prctile(dLocsSecs, 25);
        ymax = tapas_physio_prctile(dLocsSecs, 99);
        
        fs(3) = subplot(3,1,3);
        plot(t(LOCS(1:end-1)), dLocsSecs); title('duration betwenn scan events - search for bad peaks here!');
        xlabel('t (s)');
        ylabel('t (ms)');
        ylim([0.9*ymin, 1.1*ymax]);
        linkaxes(fs,'x');
        
    end
    subplot(3,1,2);
    legend(hp, lg);
    
end


%% Return error if not enough events flund
% VOLLOCS = find(abs(diff(z2))>thresh.vol);
if isempty(VOLLOCS) || isempty(LOCS)
    error('No volume start events found, Decrease thresh.vol or thresh.slice after considering the Thresholding figure');
elseif length(LOCS) < nSlicesPerBeat
    error('Too few slice start events found. Decrease thresh.slice after considering the Thresholding figure');
end

if doCountSliceEventsFromLogfileStart
    if length(VOLLOCS)< (Nprep+nScans+nDummies)
        error(['Not enough volume events found. \n\tFound:  %d\n ' ...
            '\tNeeded: %d+%d+%d (Nprep+nDummies+nScans)\n' ...
            'Please lower thresh.vol or thresh.vol_spacing\n'], ...
            length(VOLLOCS), Nprep, nDummies, nScans);
    end
else
    if length(VOLLOCS)< (nScans+nDummies)
        error(['Not enough volume events found. \n\tFound:  %d\n ' ...
            '\tNeeded: %d+%d (nDummies+nScans)\n' ...
            'Please lower thresh.vol or thresh.vol_spacing\n'], ...
            length(VOLLOCS), nDummies, nScans);
    end
end

