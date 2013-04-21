function [R, ons_secs] = main_ECG7T_manualPeakPick(files, ...
    thresh, sqpar, order, verbose)
% RETROICOR - regressor creation based on Glover, G. MRM 44, 2000
%
% USAGE
%   [R, ons_secs] = main_create_RETROICOR_regressors(logfile, ...
%    thresh, sqpar, order, verbose, regressordir, logfile_w_scanevents, orthogonalize_cardiac, rpdir)
%
%------------------------------------------------------------------------
% IN
%   files   is a structure containing the following filenames (with full
%           path)
%       .vendor             'Philips', 'GE' or 'Siemens', depending on your
%                           MR Scanner system
%       .log_cardiac        contains ECG or pulse oximeter time course
%                           for Philips: 'SCANPHYSLOG<DATE&TIME>.log';
%                           can be found on scanner in G:/log/scanphyslog-
%                           directory, one file is created per scan, make sure to take
%                           the one with the time stamp corresponding to your PAR/REC
%                           files
%       .log_respiration    contains breathing belt amplitude time course
%                           for Philips: same as .log_cardiac
%       .input_other_multiple_regressors = '';
%                           other regressors which should end up in "multiple regressors"
%                           slot of SPM-GLM; either txt-file or mat-file with variable R
%                           e.g. realignment parameters rp_*.txt
%       .output_multiple_regressors
%                           output .mat-file containing a variable R with
%                           all RETROICOR-regressors; can be inserted
%                           directly as "multiple regressors" for SPM
%                           1st level design specification
%                           e.g. 'multiple_regressors.mat' in SPM-analysis
%                           folder
%
%   thresh  - thresh is a structure with the following elements
%
%           thresh.scan_timing.
%           .zero    - gradient values below this value are set to zero;
%                      should be those which are unrelated to slice acquisition start
%           .slice   - minimum gradient amplitude to be exceeded when a slice
%                      scan starts
%           .vol     - minimum gradient amplitude to be exceeded when a new
%                      volume scan starts;
%                      leave [], if volume events shall be determined as
%                      every Nslices-th scan event
%           .grad_direction
%                    - leave empty to use nominal timing;
%                      if set, sequence timing is calculated from logged gradient timecourse;
%                    - value determines which gradient direction timecourse is used to
%                      identify scan volume/slice start events ('x', 'y', 'z')
%           .vol_spacing
%                   -  duration (in seconds) from last slice acq to
%                      first slice of next volume;
%                      leave [], if .vol-threshold shall be used
%
%           thresh.cardiac.
%
%           .modality - 'ecg' or 'oxy'; ECG or Pulse oximeter used?
%           .min -     - for modality 'ECG': [percent peak height of sample QRS wave]
%                      if set, ECG heartbeat event is calculated from ECG
%                      timeseries by detecting local maxima of
%                      cross-correlation to a sample QRS-wave;
%                      leave empty, to use Philips' log of heartbeat event
%                      - for modality 'OXY': [peak height of pulse oxymeter] if set, pulse
%                      oxymeter data is used and minimal peak height
%                      set is used to determined maxima
%           .kRpeakfile
%                    - [default: not set] string of file containing a
%                      reference ECG-peak
%                      if set, ECG peak detection via cross-correlation (via
%                      setting .ECG_min) performed with a saved reference ECG peak
%                      This file is saved after picking the QRS-wave
%                      manually (i.e. if .ECG_min is set), so that
%                      results are reproducible
%           .manual_peak_select
%                      [false] or true; if true, a user input is
%                      required to specify a characteristic R-peak interval in the ECG
%                      or pulse oximetry time series
%
%       thresh.respiratory.
%
%           .resp_max
%                   -  if set, all peaks above that breathing belt amplitude
%                      are ignored for respiratory phase evaluation
%
% NOTE: estimate gradient thresholds from visual inspection of the gradient timecourses
%       They only have to be set once per sequence, i.e. can be used for
%       all subjects and sessions
%
%   sqpar                   - sequence timing parameters
%           .Nslices        - number of slices per volume in fMRI scan
%           .NslicesPerBeat - usually equals Nslices, unless you trigger with the heart beat
%           .TR             - repetition time in seconds
%           .Ndummies       - number of dummy volumes
%           .Nscans         - number of full volumes saved (volumes in nifti file,
%                             usually rows in your design matrix)
%           .Nprep          - number of non-dummy, volume like preparation pulses
%                             before 1st dummy scan. If set, logfile is read from beginning,
%                             otherwise volumes are counted from last detected volume in the logfile
%           .TimeSliceToSlice - time between the acquisition of 2 subsequent
%                             slices; typically TR/Nslices or
%                             minTR/Nslices, if minimal temporal slice
%                             spacing was chosen
%            onset_slice    - slice whose scan onset determines the adjustment of the
%                             regressor timing to a particular slice for the whole volume
%                             NOTE: only necessary, if thresh.grad_direction is empty
%
% order     - order of RETROICOR expansion, taken from Harvey2008, JRMI28(6), p1337ff.
%       .c  - cardiac [default = 3]
%       .r  - respiratory [default = 4]
%       .cr - multiplicative terms: cardiac X respiratory [default 1]
%       .orthogonalise
%           - string indicating which regressors shall be
%             orthogonalised; mainly needed, if
%           acquisition was triggered to heartbeat (set to 'cardiac') OR
%           if session mean shall be evaluated (e.g. SFNR-studies, set to
%           'all')
%             'n' or 'none'     - no orthogonalisation is performed
%             'c' or 'cardiac'  - only cardiac regressors are orthogonalised
%             'r' or 'resp'     - only respiration regressors are orthogonalised
%             'mult'            - only multiplicative regressors are orthogonalised
%             'all'             - all physiological regressors are
%                                 orthogonalised to each other
% verbose   - create informative plots (1= yes, 0 = no)
%
%------------------------------------------------------------------------
% OUT
%   R         - matrix of physiological regressors to be plugged in as
%           "multiple regressors" into SPM design matrix
%           If realignment parameters were found in <regressordir>, they are
%           appended as columns to <R>
%   ons_secs  - onset time stamps of scan events, heartbeats, breathing (in
%           seconds)
%
% -------------------------------------------------------------------------
% Lars Kasper, August 2011
% $Id: main_create_RETROICOR_regressors.m 164 2013-02-28 11:47:07Z kasperla $
%


%% 0. set default parameters
% order of RETROICOR expansion, taken from Harvey2008, JRMI28(6), p1337ff.
if nargin < 4 || isempty(order)
    order.c = 3; %cardiac
    order.r = 4; % respiratory
    order.cr = 1; % cardiac X respiratory
    order.orthogonalise = 'none';
end

if nargin < 5, verbose = true; end;

manualCorrect = 1; % !!! include this somewhere in the input of the function later.


%% 1. Read in vendor-specific physiological log-files
[ons_secs.c, ons_secs.r, ons_secs.t, ons_secs.cpulse] = ...
    read_physlogfiles(files, thresh.cardiac.modality);

if verbose >= 1
    plot_raw_physdata(ons_secs);
end


%% 2. Create scan timing nominally or from gradient time-course
% the latter is only necessary, if no patch is used and therefore no scan event
% triggers are written into the last column of the scanphyslog-file
if isempty(thresh.scan_timing)
    [VOLLOCS, LOCS] = create_nominal_scan_timing(ons_secs.t, sqpar);
else
    [VOLLOCS, LOCS] = create_scan_timing_from_gradients_philips( ...
        files, thresh.scan_timing, sqpar, verbose);
end

[ons_secs.svolpulse, ons_secs.spulse, ons_secs.spulse_per_vol] = ...
    get_onsets_from_locs(ons_secs.t, VOLLOCS, LOCS, sqpar, verbose);


%% 3. Extract and check physiological parameters (onsets)
% plot whether physdata is alright or events are missing (too low/high
% heart rate? breathing amplitude overshoot?)

if isfield(thresh.cardiac, 'min') && ~isempty(thresh.cardiac.min)
    ons_secs.cpulse = get_cardiac_pulses(ons_secs.t, ons_secs.c, ...
        thresh.cardiac, verbose);
end

if manualCorrect
    ons_secs = physio_correct_cardiac_pulses_manually(ons_secs,80,60,50);
end

if verbose>=2
    %plot_raw_physdata_diagnostics(ons_secs.t, ons_secs.cpulse, ons_secs.r);
    [outliersHigh,outliersLow] = plot_raw_physdata_diagnostics_and_select(ons_secs.t, ons_secs.cpulse, ons_secs.r);
end

[ons_secs, sqpar] = crop_scanphysevents_to_acq_window(ons_secs, sqpar);
if verbose >= 1
    plot_cropped_phys_to_acqwindow(ons_secs, sqpar);
end


%% 4. Create RETROICOR regressors for SPM
[cardiac_sess, respire_sess, mult_sess] = ...
    create_retroicor_regressors(ons_secs, sqpar, thresh, ...
    order, verbose);


% 4.1.  Load other confound regressors, e.g. realigment parameters
if isfield(files, 'input_other_multiple_regressors') && ~isempty(files.input_other_multiple_regressors)
    input_R = load_other_multiple_regressors(files.input_other_multiple_regressors);
else
    input_R = [];
end


% 4.2   Orthogonalisation of regressors ensures numerical stability for
%       otherwise correlated cardiac regressors
R = orthogonalise_physiological_regressors(cardiac_sess, respire_sess, ...
    mult_sess, input_R, order.orthogonalise, verbose);

% 4.3   Save Multiple Regressors file for SPM
% TODO: save as TXT-file for non-Matlab use
save(files.output_multiple_regressors, 'R');