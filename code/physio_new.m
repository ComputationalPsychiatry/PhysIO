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
%   phzsio         - the complete physio structure, which can be unsed in
%                     physio_main_create_regressors
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
%
%
%
% model
% .type  'RETROICOR'
% .input_other_multiple_regressors = '';
%                           other regressors which should end up in "multiple regressors"
%                           slot of SPM-GLM; either txt-file or mat-file with variable R
%                           e.g. realignment parameters rp_*.txt
% .output_multiple_regressors
%                           output .mat-file containing a variable R with
%                           all RETROICOR-regressors; can be inserted
%                           directly as "multiple regressors" for SPM
%                           1st level design specification
%                           e.g. 'multiple_regressors.mat' in SPM-analysis
%                           folder
% .order     - order of RETROICOR expansion, taken from Harvey2008, JRMI28(6), p1337ff.
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
%           .Nprep          - number of non-dummy, volume like preparation pulses
%                             before 1st dummy scan. If set, logfile is read from beginning,
%                             otherwise volumes are counted from last detected volume in the logfile
sqpar.TimeSliceToSlice  = [];
%           .TimeSliceToSlice - time between the acquisition of 2 subsequent
%                             slices; typically TR/Nslices or
%                             minTR/Nslices, if minimal temporal slice
%                             spacing was chosen
%                             NOTE: only necessary, if
%                             thresh.grad_direction is empty and nominal
%                             scan timing is used
    sqpar.onset_slice       = 19;
%            .onset_slice    - slice whose scan onset determines the adjustment of the
%                             regressor timing to a particular slice for the whole volume
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
    
    thresh.cardiac.initial_cpulse_select.method = 'load_from_logfile'; % 'load_from_logfile', 'manual', 'load'
    thresh.cardiac.initial_cpulse_select.file = '';
    thresh.cardiac.initial_cpulse_select.min = [];    
    thresh.cardiac.initial_cpulse_select.kRpeak = [];
    
    thresh.cardiac.posthoc_cpulse_select.method = 'off'; % 'off', 'manual', 'load'
    thresh.cardiac.posthoc_cpulse_select.file = '';        
    thresh.cardiac.posthoc_cpulse_select.percentile = 80;
    thresh.cardiac.posthoc_cpulse_select.upperThresh = 60;
    thresh.cardiac.posthoc_cpulse_select.lowerThresh = 60; 
    
    %% order
    model.type = '';
    model.input_other_multiple_regressors = ''; % either txt-file or mat-file with variable R
    model.output_multiple_regressors = '';
    model.order = struct('c',[],'r',[],'cr',[], 'orthogonalise', '');

    %% verbose
    verbose.level = [];
    verbose.fig_handles = [];
    verbose.fig_output_file = '';
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
