function physio = physio_new(default_scheme, physio_in)
% creates complete PhysIO structure to be fed into physio_main_create_regressors
%   output = physio_new(input)
%
% IN
%   default_scheme - if set, default values for structure entries are set
%                       according to the application
%                       different templates are predefined, e.g.
%                       'empty' (default) - all strings are set to '', all
%                                     numbers to []
%                       'RETROICOR'
%                   `       % order of RETROICOR expansion taken from Harvey2008, JRMI28(6), p1337ff.
%                       'scan_timing_from_start'
%                       'manual_peak_select'
%   physio_in       - used as input, only the fields related to the default_scheme
%                     are overwritten, the others are kept as in physio_in
% OUT
%   physio         - the complete physio structure, which can be unsed in
%                     physio_main_create_regressors
% NOTE
%   All parameters used in the physIO toolbox are defined AND DOCUMENTED in
%   this file. Just scroll down and read through the comments!
% 
% EXAMPLE
%   physio = physio_new('empty')
%   physio = physio_new('RETROICOR');
%   physio = physio_new('manual_peak_select', physio);
%
%   See also physio_main_create_regressors
%
% Author: Lars Kasper
% Created: 2013-04-23
% Copyright (C) 2013 TNU, Institute for Biomedical Engineering, University of Zurich and ETH Zurich.
%
% This file is part of the TNU CheckPhysRETROICOR toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id$
%
%
% if not specified differently, create everything empty
if ~nargin
    default_scheme = 'empty';
end

if nargin >= 2
    log_files   = physio_in.log_files;
    thresh  = physio_in.thresh;
    sqpar   = physio_in.sqpar;
    model   = physio_in.model;
    verbose = physio_in.verbose;
else
    
    
    %% log_files
    % structure containing general physiological log-file information
    log_files.vendor       = ''; % 'Philips', 'GE', or 'Siemens', depending on your
                                 %  MR Scanner system
    log_files.cardiac      = ''; % 'SCANPHYSLOG.log'; logfile with cardiac data
    log_files.respiration  = ''; % 'SCANPHYSLOG.log'; logfile with respiratory data
                                 %                    (same as .cardiac for Philips)
  % log_files.sampling_interval = [];   % in seconds, 2e-3 for Philips, variable for GE,
                                        % e.g. 40e-3
 
                                        
    %% sqpar
    sqpar.Nslices           = [];   % number of slices per volume in fMRI scan
    sqpar.NslicesPerBeat    = [];   % usually equals Nslices, unless you trigger with the heart beat
    sqpar.TR                = [];   % volume repetition time in seconds
    sqpar.Ndummies          = [];   % number of dummy volumes
    sqpar.Nscans            = [];   % number of full volumes saved (volumes in nifti file,
                                    % usually rows in your design matrix)
    sqpar.Nprep             = [];   % set to >=0 to count scans and dummy
                                    % number of non-dummy, volume like preparation pulses
                                    % before 1st dummy scan. If set, logfile is read from beginning,
                                    % otherwise volumes are counted from last detected volume in the logfile
    sqpar.TimeSliceToSlice  = [];   % time between the acquisition of 2 subsequent
                                    % slices; typically TR/Nslices or minTR/Nslices, 
                                    % if minimal temporal slice spacing was chosen
                                    % NOTE: only necessary, if thresh.grad_direction 
                                    % is empty and nominal scan timing is used
    sqpar.onset_slice       = 19;   % slice whose scan onset determines the adjustment of the
                                    % regressor timing to a particular slice for the whole volume
                                    % volumes from beginning of run, i.e. logfile,
                                    % includes counting of preparation gradients
                              
    
    %% thresh
    % determines thresholds used in preprocessing physiological logfiles,
    % either their timing (thresh.scan_timing) or the peripheral measures
    % itself (thresh.cardiac, thresh.respiration)
    thresh.scan_timing = [];    % leave empty, if nominal scan timing, 
                                % derived from sqpar, shall be used

    thresh.scan_timing.grad_direction = ''; % 'x', 'y', or 'z'; 
                                            % if set, sequence timing is calculated 
                                            % from logged gradient timecourse along 
                                            % this coordinate axis;
    thresh.scan_timing.zero     = [];   % gradient values below this value are set to zero;
                                        % should be those which are unrelated to slice acquisition start
    thresh.scan_timing.slice    = [];   % minimum gradient amplitude to be exceeded when a slice scan starts
    thresh.scan_timing.vol      = [];   % minimum gradient amplitude to be exceeded when a new
                                        % volume scan starts; 
                                        % leave [], if volume events shall be determined as
                                        % every Nslices-th scan event or via vol_spacing
    thresh.vol_spacing          = [];   % duration (in seconds) from last slice acq to
                                        % first slice of next volume;
                                        % leave [], if .vol-threshold shall be used
                                        
    thresh.cardiac = [];
    thresh.cardiac.modality = ''; % 'ECG','ECG_raw', or 'OXY' (for pulse oximetry), 'OXY_OLD', [deprecated]
    
    % The initial cardiac pulse selection structure: Determines how the
    % majority of cardiac pulses is detected
    thresh.cardiac.initial_cpulse_select.method = 'load_from_logfile'; % 'load_from_logfile', 'manual' (rather: threshold...autocorrelate?), 'load'
    thresh.cardiac.initial_cpulse_select.file = ''; % file containing reference ECG-peak (variable kRpeak)
                                                    % used for method 'manual' or 'load' [default: not set] string of file containing a
                                                    % if method == 'manual', this file is saved after picking the QRS-wave
                                                    % such that results are reproducible                                    
    thresh.cardiac.initial_cpulse_select.min = [];  % threshold for correlation with QRS-wave to find cardiac pulses 
    thresh.cardiac.initial_cpulse_select.kRpeak = []; % variable saving an example cardiac QRS-wave to correlate with ECG time series
    
    % The posthoc cardiac pulse selection structure: If only few (<20)
    % cardiac pulses are missing in a session due to bad signal quality, a
    % manual selection after visual inspection is possible using the
    % following parameters. The results are saved for reproducibility
    thresh.cardiac.posthoc_cpulse_select.method = 'off'; % 'off', 'manual', 'load'
                                                         % 'off' - no manual selection of peaks
                                                         % 'manual' - pick and save additional peaks manually
                                                         % 'load' - load previously selected cardiac pulses                                                          
    thresh.cardiac.posthoc_cpulse_select.file = '';  % filename where cardiac pulses are saved after manual picking      
    
    % Suspicious positions of missing or too many cardiac pulses are
    % pre-selected by detecting outliers in histogram of
    % heart-beat-2-beat-intervals
    thresh.cardiac.posthoc_cpulse_select.percentile = 80; % percentile of beat-2-beat interval histogram that constitutes the "average heart beat duration" in the session
    thresh.cardiac.posthoc_cpulse_select.upperThresh = 60; % minimum exceedance (in %) from average heartbeat duration to be classified as missing heartbeat
    thresh.cardiac.posthoc_cpulse_select.lowerThresh = 60; % minimum reduction (in %) from average heartbeat duration to be classified an abundant heartbeat
    
    
    %% model
    % Determines the physiological noise model derived from preprocessed physiological data
    model.type = '';                            % 'RETROICOR' - as in Glover el al, MRM 44, 2000
    model.input_other_multiple_regressors = ''; % other nuisance regressors to be included in design matrix
                                                % either txt-file or mat-file with variable R
    model.output_multiple_regressors = '';      % output file for usage in SPM multiple_regressors GLM-specification
                                                % either txt-file or mat-file with variable R
    model.order.c = [];                         % natural number, order of cardiac phase Fourier expansion
    model.order.r = [];                         % natural number, order of respiratory phase Fourier expansion
    model.order.cr = [];                        % natural number, order of cardiac-respiratory-phase-interaction Fourier expansion
                                                % See Harvey et al, JMRI 28, 2008
    model.order.orthogonalise = 'none';         % string indicating which regressors shall be orthogonalised; 
                                                % mainly needed, if acquisition was triggered to heartbeat (set to 'cardiac') OR
                                                % if session mean shall be evaluated (e.g. SFNR-studies, set to 'all')
                                                % 'n' or 'none'     - no orthogonalisation is performed
                                                % Possible Values (default: 'none'
                                                %   'c' or 'cardiac'  - only cardiac regressors are orthogonalised
                                                %   'r' or 'resp'     - only respiration regressors are orthogonalised
                                                %   'mult'            - only multiplicative regressors are orthogonalised
                                                %   'all'             - all physiological regressors are orthogonalised to each other

                                                
    %% verbose
    % determines how many figures shall be generated to follow the workflow
    % of the toolbox and whether the graphical output shall be saved (to a
    % PostScript-file)
    verbose.level = 1;            % 0 = no graphical output; 1 = main plots (default);  
                                  % 2 = debugging plots, for setting up new study; 3 = all plots
    verbose.fig_handles = [];     % collector of all generated figure handles during a run of physio_main_create_regressors
    verbose.fig_output_file = ''; % file name (including extension) where to print all physIO output figures to,
                                  % e.g. 'PhysIO_output.ps' or 'PhysIO_output.jpg'
                                  % The specified extension determines how the
                                  % figures will be saved
                                  %     .ps - all figures are saved to the
                                  %     same, multiple-page postscript-file
                                  %     .fig, .tiff,  .jpg 
                                  %         - one file is created for each
                                  %         figure, appended by its figure
                                  %         index, e.g. 'PhysIO_output_fig01.jpg'
end

switch default_scheme
    case 'RETROICOR'
        model.type = 'RETROICOR';
        model.order = struct('c',3,'r',4,'cr',1, 'orthogonalise', 'none');
    case 'redetect_peaks_from_logfile'
        thresh.cardiac.initial_cpulse_select.method = 'manual'; % 'load_from_logfile', 'manual', 'load'
        thresh.cardiac.initial_cpulse_select.file = 'kRpeak.mat';
        thresh.cardiac.initial_cpulse_select.min = 1;
        thresh.cardiac.initial_cpulse_select.kRpeak = [];
    case 'manual_peak_select'
        thresh.cardiac.posthoc_cpulse_select.method = 'manual'; % 'off', 'manual' or 'load',
        thresh.cardiac.posthoc_cpulse_select.file = 'posthoc_cpulse.mat';
        thresh.cardiac.posthoc_cpulse_select.percentile = 80;
        thresh.cardiac.posthoc_cpulse_select.upperThresh = 60;
        thresh.cardiac.posthoc_cpulse_select.lowerThresh = 30;        
end

%% assemble output
physio.log_files   = log_files;
physio.thresh  = thresh;
physio.sqpar   = sqpar;
physio.model   = model;
physio.verbose = verbose;
