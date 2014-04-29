function physio = tapas_physio_cfg_matlabbatch
% Lars Kasper, March 2013
%
% Copyright (C) 2013, Institute for Biomedical Engineering, ETH/Uni Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id$


pathThis = fileparts(mfilename('fullpath')); % TODO: more elegant via SPM!
addpath(pathThis);

%==========================================================================
%% Sub-structure log_files
%==========================================================================

%--------------------------------------------------------------------------
% vendor
%--------------------------------------------------------------------------
vendor        = cfg_menu;
vendor.tag    = 'vendor';
vendor.name   = 'vendor';
vendor.help   = {'Choose Vendor of your scanner from list or Custom'
    'Custom logfiles should be ASCII-files with one sample per row'};
vendor.labels = {'Philips', 'GE', 'Siemens', 'Custom'};
vendor.values = {'Philips', 'GE', 'Siemens', 'Custom'};
vendor.val    = {'Philips'};

%--------------------------------------------------------------------------
% cardiac
%--------------------------------------------------------------------------
cardiac         = cfg_files;
cardiac.tag     = 'cardiac';
cardiac.name    = 'log_cardiac';
%cardiac.val     = {{'/Users/kasperla/Documents/code/matlab/smoothing_trunk/tSNR_fMRI_SPM/CheckPhysRETROICOR/PhysIOToolbox/examples/Philips/ECG3T/SCANPHYSLOG.log'}};
cardiac.help    = {'...'};
cardiac.filter  = 'any';
cardiac.ufilter = '.*';
cardiac.num     = [1 1];

%--------------------------------------------------------------------------
% respiration (filename)
%--------------------------------------------------------------------------
respiration         = cfg_files;
respiration.tag     = 'respiration';
respiration.name    = 'log_respiration';
% respiration.val     = {{'/Users/kasperla/Documents/code/matlab/smoothing_trunk/tSNR_fMRI_SPM/CheckPhysRETROICOR/PhysIOToolbox/examples/Philips/ECG3T/SCANPHYSLOG.log'}};
respiration.help    = {'...'};
respiration.filter  = 'any';
respiration.ufilter = '.*';
respiration.num     = [0 1];


%--------------------------------------------------------------------------
% sampling_interval
%--------------------------------------------------------------------------
sampling_interval         = cfg_entry;
sampling_interval.tag     = 'sampling_interval';
sampling_interval.name    = 'sampling_interval';
sampling_interval.help    = {
    'sampling interval of phys log files (in seconds)'
    ' If empty, default values are used: 2 ms for Philips, 25 ms for GE and others'
};
sampling_interval.strtype = 'e';
sampling_interval.num     = [Inf Inf];
sampling_interval.val     = {[]};

%--------------------------------------------------------------------------
% relative_start_acquisition
%--------------------------------------------------------------------------
relative_start_acquisition         = cfg_entry;
relative_start_acquisition.tag     = 'relative_start_acquisition';
relative_start_acquisition.name    = 'relative_start_acquisition';
relative_start_acquisition.help    = {'start time of 1st scan (or dummy) relative to start of phys logfile'};
relative_start_acquisition.strtype = 'e';
relative_start_acquisition.num     = [Inf Inf];
relative_start_acquisition.val     = {0};


%--------------------------------------------------------------------------
% files
%--------------------------------------------------------------------------
files      = cfg_branch;
files.tag  = 'log_files';
files.name = 'log_files';
files.val  = {vendor cardiac respiration sampling_interval, relative_start_acquisition};
files.help = {'...'};




%==========================================================================
%% Sub-structure sqpar
%==========================================================================



%--------------------------------------------------------------------------
% Nscans
%--------------------------------------------------------------------------
Nscans         = cfg_entry;
Nscans.tag     = 'Nscans';
Nscans.name    = 'Nscans';
Nscans.help    = {'Number of scans (volumes) in design matrix'};
Nscans.strtype = 'e';
Nscans.num     = [Inf Inf];
%Nscans.val     = {495};

%--------------------------------------------------------------------------
% Ndummies
%--------------------------------------------------------------------------
Ndummies         = cfg_entry;
Ndummies.tag     = 'Ndummies';
Ndummies.name    = 'Ndummies';
Ndummies.help    = {
    'Number of dummies that were acquired (but do not show up in design matrix'
    '(also enter correct number, if dummies are not saved in imaging file'
    };
Ndummies.strtype = 'e';
Ndummies.num     = [Inf Inf];
%Ndummies.val     = {3};

%--------------------------------------------------------------------------
% TR
%--------------------------------------------------------------------------
TR         = cfg_entry;
TR.tag     = 'TR';
TR.name    = 'TR';
TR.help    = {'Repetition time in seconds'};
TR.strtype = 'e';
TR.num     = [Inf Inf];
%TR.val     = {2.5};

%--------------------------------------------------------------------------
% NslicesPerBeat
%--------------------------------------------------------------------------
NslicesPerBeat         = cfg_entry;
NslicesPerBeat.tag     = 'NslicesPerBeat';
NslicesPerBeat.name    = 'NslicesPerBeat';
NslicesPerBeat.help    = {'Only for triggered (gated) sequences: '
    'Number of slices acquired per heartbeat'};
NslicesPerBeat.strtype = 'e';
NslicesPerBeat.num     = [Inf Inf];
NslicesPerBeat.val     = {[]};


%--------------------------------------------------------------------------
% Nslices
%--------------------------------------------------------------------------
Nslices         = cfg_entry;
Nslices.tag     = 'Nslices';
Nslices.name    = 'Nslices';
Nslices.help    = {'Number of slices in one volume'};
Nslices.strtype = 'e';
Nslices.num     = [Inf Inf];
%Nslices.val     = {37};



%--------------------------------------------------------------------------
% onset_slice
%--------------------------------------------------------------------------
onset_slice         = cfg_entry;
onset_slice.tag     = 'onset_slice';
onset_slice.name    = 'onset_slice';
onset_slice.help    = {
    'slice to which regressors are temporally aligned'
    'Typically the slice where your most important activation is expected'};
onset_slice.strtype = 'e';
onset_slice.num     = [Inf Inf];
%onset_slice.val     = {19};

%--------------------------------------------------------------------------
% Nprep
%--------------------------------------------------------------------------
Nprep         = cfg_entry;
Nprep.tag     = 'Nprep';
Nprep.name    = 'Nprep';
Nprep.help    = {'Preparation (e.g. shimming) volumes acquired before first dummy'};
Nprep.strtype = 'e';
Nprep.num     = [Inf Inf];
Nprep.val     = {[]};

%--------------------------------------------------------------------------
% time_slice_to_slice
%--------------------------------------------------------------------------
time_slice_to_slice         = cfg_entry;
time_slice_to_slice.tag     = 'time_slice_to_slice';
time_slice_to_slice.name    = 'time_slice_to_slice';
time_slice_to_slice.help    = {
    'duration between acquisition of two different slices'
    'if empty, set to default value TR/Nslices'
    'differs e.g. if slice timing was minimal and TR was bigger than needed'
    'to acquire Nslices'
    };
time_slice_to_slice.strtype = 'e';
time_slice_to_slice.num     = [Inf Inf];
time_slice_to_slice.val     = {[]};

%--------------------------------------------------------------------------
% sqpar
%--------------------------------------------------------------------------
sqpar      = cfg_branch;
sqpar.tag  = 'sqpar';
sqpar.name = 'sqpar (Sequence timing parameters)';
sqpar.val  = {Nslices NslicesPerBeat TR Ndummies Nscans onset_slice time_slice_to_slice Nprep};
sqpar.help = {'...'};




%==========================================================================
%% Sub-structure model
%==========================================================================


%--------------------------------------------------------------------------
% c
%--------------------------------------------------------------------------
c         = cfg_entry;
c.tag     = 'c';
c.name    = 'cardiac';
c.help    = {'...'};
c.strtype = 'e';
c.num     = [1 1];
c.val     = {3};

%--------------------------------------------------------------------------
% r
%--------------------------------------------------------------------------
r         = cfg_entry;
r.tag     = 'r';
r.name    = 'respiratory';
r.help    = {'...'};
r.strtype = 'e';
r.num     = [1 1];
r.val     = {4};

%--------------------------------------------------------------------------
% cr
%--------------------------------------------------------------------------
cr         = cfg_entry;
cr.tag     = 'cr';
cr.name    = 'cardiac X respiratory';
cr.help    = {'...'};
cr.strtype = 'e';
cr.num     = [1 1];
cr.val     = {1};

%--------------------------------------------------------------------------
% orthog
%--------------------------------------------------------------------------
orthog        = cfg_menu;
orthog.tag    = 'orthogonalise';
orthog.name   = 'orthogonalise';
orthog.help   = {'...'};
orthog.labels = {'none' 'cardiac' 'resp' 'mult' 'all'};
orthog.values = {'none' 'cardiac' 'resp' 'mult' 'all'};
orthog.val    = {'none'};

%--------------------------------------------------------------------------
% order
%--------------------------------------------------------------------------
order      = cfg_branch;
order.tag  = 'order';
order.name = 'order';
order.val  = {c r cr orthog};
order.help = {'...'};

%--------------------------------------------------------------------------
% model_type
%--------------------------------------------------------------------------
model_type        = cfg_menu;
model_type.tag    = 'type';
model_type.name   = 'type';
model_type.help   = {'...'};
model_type.labels = {
    'RETROICOR (RETRO)' 
    'Heart Rate Variability (HRV)'
    'Respiratory Volume per Time (RVT)'
    'RETRO+HRV'
    'RETRO+RVT'
    'HRV+RVT'
    'RETRO+HRV+RVT'
    };
model_type.values = {
    'RETROICOR' 
    'HRV'
    'RVT'
    'RETROICOR_HRV'
    'RETROICOR_RVT'
    'HRV_RVT'
    'RETROICOR_HRV_RVT'
    };
model_type.val    = {'RETROICOR'};

%--------------------------------------------------------------------------
% output_multiple_regressors
%--------------------------------------------------------------------------
output_multiple_regressors         = cfg_entry;
output_multiple_regressors.tag     = 'output_multiple_regressors';
output_multiple_regressors.name    = 'output_multiple_regressors';
output_multiple_regressors.help    = {'...'};
output_multiple_regressors.strtype = 's';
output_multiple_regressors.num     = [1 Inf];
output_multiple_regressors.val     = {'multiple_regressors.txt'};

%--------------------------------------------------------------------------
% input_other_multiple_regressors
%--------------------------------------------------------------------------
input_other_multiple_regressors         = cfg_files;
input_other_multiple_regressors.tag     = 'input_other_multiple_regressors';
input_other_multiple_regressors.name    = 'input_other_multiple_regressors';
input_other_multiple_regressors.val     = {{''}};
input_other_multiple_regressors.help    = {'...'};
input_other_multiple_regressors.filter  = '.*';
input_other_multiple_regressors.ufilter = '.mat$|.txt$';
input_other_multiple_regressors.num     = [0 1];

%--------------------------------------------------------------------------
% model
%--------------------------------------------------------------------------
model      = cfg_branch;
model.tag  = 'model';
model.name = 'model';
model.val  = {model_type, order, input_other_multiple_regressors, ...
    output_multiple_regressors};
model.help = {'...'};




% ==========================================================================
%% Sub-structure thresh
%==========================================================================


% ==========================================================================
%% Subsub-structure scan_timing
%==========================================================================

%--------------------------------------------------------------------------
% scan_timing_method
%--------------------------------------------------------------------------
scan_timing_method        = cfg_menu;
scan_timing_method.tag    = 'method';
scan_timing_method.name   = 'method';
scan_timing_method.help   = {
 'method to determine slice onset times for regressors'
'''nominal'' - to derive slice acquisition timing from sqpar directly'
'''gradient'' or ''gradient_log'' - derive from logged gradient time courses'
    };
scan_timing_method.labels = {'nominal' 'gradient_log'};
scan_timing_method.values = {'nominal' 'gradient_log'};
scan_timing_method.val    = {'nominal'};


%--------------------------------------------------------------------------
% grad_direction
%--------------------------------------------------------------------------
grad_direction        = cfg_menu;
grad_direction.tag    = 'grad_direction';
grad_direction.name   = 'grad_direction';
grad_direction.help   = {'...'};
grad_direction.labels = {'x' 'y' 'z'};
grad_direction.values = {'x' 'y' 'z'};
grad_direction.val    = {'y'};

%--------------------------------------------------------------------------
% vol_spacing
%--------------------------------------------------------------------------
vol_spacing         = cfg_entry;
vol_spacing.tag     = 'vol_spacing';
vol_spacing.name    = 'vol_spacing';
vol_spacing.help    = {'time (in ms) between last slice of n-th volume' 
    'and 1st slice of n+1-th volume(overrides .vol-threshold)'
    'NOTE: Leave empty if .vol shall be used'};
vol_spacing.strtype = 'e';
vol_spacing.num     = [Inf Inf];
vol_spacing.val     = {[]};

%--------------------------------------------------------------------------
% vol
%--------------------------------------------------------------------------
vol         = cfg_entry;
vol.tag     = 'vol';
vol.name    = 'vol';
vol.help    = {'Gradient Amplitude Threshold for Start of new Volume'};
vol.strtype = 'e';
vol.num     = [Inf Inf];
vol.val     = {[]};

%--------------------------------------------------------------------------
% slice
%--------------------------------------------------------------------------
slice         = cfg_entry;
slice.tag     = 'slice';
slice.name    = 'slice';
slice.help    = {'Gradient Amplitude Threshold for Start of new slice'};
slice.strtype = 'e';
slice.num     = [Inf Inf];
slice.val     = {1800};

%--------------------------------------------------------------------------
% zero
%--------------------------------------------------------------------------
zero         = cfg_entry;
zero.tag     = 'zero';
zero.name    = 'zero';
zero.help    = {'Gradient Amplitude Threshold below which values will be set to 0.'};
zero.strtype = 'e';
zero.num     = [Inf Inf];
zero.val     = {1700};

%--------------------------------------------------------------------------
% scan_timing
%--------------------------------------------------------------------------
scan_timing      = cfg_branch;
scan_timing.tag  = 'scan_timing';
scan_timing.name = 'scan_timing';
scan_timing.val  = {scan_timing_method grad_direction zero slice vol vol_spacing};
scan_timing.help = {'Determines scan timing from nominal scan parameters or logged gradient time courses'};


% ==========================================================================
%% Subsub-structure cardiac
%==========================================================================


%--------------------------------------------------------------------------
% modality
%--------------------------------------------------------------------------
modality        = cfg_menu;
modality.tag    = 'modality';
modality.name   = 'modality';
modality.help   = {'Shall ECG or PPU data be read from logfiles?'};
modality.labels = {'ECG', 'OXY/PPU'};
modality.values = {'ECG', 'PPU'};
modality.val    = {'ECG'};


%--------------------------------------------------------------------------
% initial_cpulse_select_method
%--------------------------------------------------------------------------
initial_cpulse_select_method        = cfg_menu;
initial_cpulse_select_method.tag    = 'method';
initial_cpulse_select_method.name   = 'method';
initial_cpulse_select_method.help   = {
     'The initial cardiac pulse selection structure: Determines how the'
    'majority of cardiac pulses is detected'
    ' ''auto''    - auto generation of representative QRS-wave; detection via'
    '             maximising auto-correlation with it'
    ' ''load_from_logfile'' - from phys logfile, detected R-peaks of scanner' 
    ' ''manual''  - via manually selected QRS-wave for autocoreelations'
    ' ''load''    - from previous manual/auto run'
   };
initial_cpulse_select_method.labels = {
    'auto', 'load_from_logfile', 'manual', 'load'};
initial_cpulse_select_method.values = { 'auto', 'load_from_logfile', 'manual', 'load'};
initial_cpulse_select_method.val    = {'load_from_logfile'};

%--------------------------------------------------------------------------
% initial_cpulse_select_file
%--------------------------------------------------------------------------
initial_cpulse_select_file         = cfg_entry;
initial_cpulse_select_file.tag     = 'file';
initial_cpulse_select_file.name    = 'file';
initial_cpulse_select_file.help    = {'...'};
initial_cpulse_select_file.strtype = 's';
initial_cpulse_select_file.num     = [1 Inf];
initial_cpulse_select_file.val     = {'kRpeakfile.mat'};



%--------------------------------------------------------------------------
% min
%--------------------------------------------------------------------------
min       = cfg_entry;
min.tag     = 'min';
min.name    = 'min';
min.help    = {'minimum correlation value considered a peak (for auto, manual, load-methods).'};
min.strtype = 'e';
min.num     = [Inf Inf];
min.val     = {0.4};


%--------------------------------------------------------------------------
% initial_cpulse_select
%--------------------------------------------------------------------------
initial_cpulse_select      = cfg_branch;
initial_cpulse_select.tag  = 'initial_cpulse_select';
initial_cpulse_select.name = 'initial_cpulse_select';
initial_cpulse_select.val  = {initial_cpulse_select_method min initial_cpulse_select_file};
initial_cpulse_select.help = {
    'The initial cardiac pulse selection structure: Determines how the'
    'majority of cardiac pulses is detected.'
   };



%--------------------------------------------------------------------------
% posthoc_cpulse_select_method
%--------------------------------------------------------------------------
posthoc_cpulse_select_method        = cfg_menu;
posthoc_cpulse_select_method.tag    = 'method';
posthoc_cpulse_select_method.name   = 'method';
posthoc_cpulse_select_method.help   = {
    '''off'' - no manual selection of peaks'
    '''manual'' - pick and save additional peaks manually'
    '''load'' - load previously selected cardiac pulses'
   };
posthoc_cpulse_select_method.labels = {
    'off', 'manual', 'load'};
posthoc_cpulse_select_method.values = {'off', 'manual', 'load'};
posthoc_cpulse_select_method.val    = {'off'};

%--------------------------------------------------------------------------
% posthoc_cpulse_select_file
%--------------------------------------------------------------------------
posthoc_cpulse_select_file         = cfg_entry;
posthoc_cpulse_select_file.tag     = 'file';
posthoc_cpulse_select_file.name    = 'file';
posthoc_cpulse_select_file.help    = {'...'};
posthoc_cpulse_select_file.strtype = 's';
posthoc_cpulse_select_file.num     = [1 Inf];
posthoc_cpulse_select_file.val     = {'cpulse.mat'};



%--------------------------------------------------------------------------
% posthoc_cpulse_select_percentile
%--------------------------------------------------------------------------
posthoc_cpulse_select_percentile       = cfg_entry;
posthoc_cpulse_select_percentile.tag     = 'percentile';
posthoc_cpulse_select_percentile.name    = 'percentile';
posthoc_cpulse_select_percentile.help    = {
    'percentile of beat-2-beat interval histogram that constitutes the'
    'average heart beat duration in the session'};
posthoc_cpulse_select_percentile.strtype = 'e';
posthoc_cpulse_select_percentile.num     = [Inf Inf];
posthoc_cpulse_select_percentile.val     = {80};

%--------------------------------------------------------------------------
% posthoc_cpulse_select_upper_thresh
%--------------------------------------------------------------------------
posthoc_cpulse_select_upper_thresh       = cfg_entry;
posthoc_cpulse_select_upper_thresh.tag     = 'upper_thresh';
posthoc_cpulse_select_upper_thresh.name    = 'upper_thresh';
posthoc_cpulse_select_upper_thresh.help    = {
    'minimum exceedance (in %) from average heartbeat duration '
    'to be classified as missing heartbeat'};
posthoc_cpulse_select_upper_thresh.strtype = 'e';
posthoc_cpulse_select_upper_thresh.num     = [Inf Inf];
posthoc_cpulse_select_upper_thresh.val     = {60};

%--------------------------------------------------------------------------
% posthoc_cpulse_select_lower_thresh
%--------------------------------------------------------------------------
posthoc_cpulse_select_lower_thresh       = cfg_entry;
posthoc_cpulse_select_lower_thresh.tag     = 'lower_thresh';
posthoc_cpulse_select_lower_thresh.name    = 'lower_thresh';
posthoc_cpulse_select_lower_thresh.help    = {
    'minimum reduction (in %) from average heartbeat duration'
    'to be classified an abundant heartbeat'};
posthoc_cpulse_select_lower_thresh.strtype = 'e';
posthoc_cpulse_select_lower_thresh.num     = [Inf Inf];
posthoc_cpulse_select_lower_thresh.val     = {60};


%--------------------------------------------------------------------------
% posthoc_cpulse_select
%--------------------------------------------------------------------------
posthoc_cpulse_select      = cfg_branch;
posthoc_cpulse_select.tag  = 'posthoc_cpulse_select';
posthoc_cpulse_select.name = 'posthoc_cpulse_select';
posthoc_cpulse_select.val  = {posthoc_cpulse_select_method ...
    posthoc_cpulse_select_file ...
    posthoc_cpulse_select_percentile ...
    posthoc_cpulse_select_upper_thresh ...
    posthoc_cpulse_select_lower_thresh};
posthoc_cpulse_select.help = {
    'The posthoc cardiac pulse selection structure: If only few (<20)'
    'cardiac pulses are missing in a session due to bad signal quality, a'
    'manual selection after visual inspection is possible using the'
    'following parameters. The results are saved for reproducibility.' 
   };



%--------------------------------------------------------------------------
% cardiac
%--------------------------------------------------------------------------
cardiac      = cfg_branch;
cardiac.tag  = 'cardiac';
cardiac.name = 'cardiac';
cardiac.val  = {modality initial_cpulse_select posthoc_cpulse_select};
cardiac.help = {'...'};


%--------------------------------------------------------------------------
% thresh
%--------------------------------------------------------------------------
thresh      = cfg_branch;
thresh.tag  = 'thresh';
thresh.name = 'thresh (Thresholding parameters for de-noising and timing)';
thresh.val  = {scan_timing cardiac};
thresh.help = {'Thresholding parameters for de-noising of raw peripheral data'
    'and determination of sequence timing from logged MR gradient time courses'};



%==========================================================================
%% Sub-structure verbose
%==========================================================================

%--------------------------------------------------------------------------
% level
%--------------------------------------------------------------------------
level         = cfg_entry;
level.tag     = 'level';
level.name    = 'level';
level.help    = {'...'};
level.strtype = 'e';
level.num     = [Inf Inf];
level.val     = {2};

%--------------------------------------------------------------------------
% fig_output_file
%--------------------------------------------------------------------------
fig_output_file         = cfg_entry;
fig_output_file.tag     = 'fig_output_file';
fig_output_file.name    = 'fig_output_file';
fig_output_file.help    = {'file name where figures are saved to; leave empty to not save'};
fig_output_file.strtype = 's';
fig_output_file.num     = [1 Inf];
fig_output_file.val     = {'PhysIO_output_level2.fig'};




%--------------------------------------------------------------------------
% use_tabs
%--------------------------------------------------------------------------
use_tabs        = cfg_menu;
use_tabs.tag    = 'use_tabs';
use_tabs.name   = 'use_tabs';
use_tabs.help   = {'use spm_tabs for plotting'};
use_tabs.labels = {'true' 'false'};
use_tabs.values = {true, false};
use_tabs.val    = {false};


%--------------------------------------------------------------------------
% verbose
%--------------------------------------------------------------------------
verbose        = cfg_branch;
verbose.tag    = 'verbose';
verbose.name   = 'verbose';
verbose.help   = {
' determines how many figures shall be generated to follow the workflow'
    ' of the toolbox and whether the graphical output shall be saved (to a'
    ' PostScript-file)'
    ' 0 = no graphical output;'
    ' 1 = (default) main plots : Fig 1: gradient scan timing (if selected) ;'
    '                            Fig 2: heart beat/breathing statistics & outlier;'
    '                            Fig 3: final multiple_regressors matrix'
    ' 2 = debugging plots        for setting up new study or if Fig 2 had'
    '                            outliers'
    '                            Fig 1: raw phys logfile data'
    '                            Fig 2: gradient scan timing (if selected)'
    '                            Fig 3: cutout interval of logfile for'
    '                            regressor creation (including scan timing'
    '                            and raw phys data)'
    '                            Fig 4: heart beat/breathing statistics & outlier;'
    '                            Fig 5: time course of all sampled RETROICOR'
    '                                   regressors'
    '                            Fig 6: final multiple_regressors matrix'
    ''
    ' 3 = all plots'
    '                            Fig 1: raw phys logfile data'
    '                            Fig 2: gradient scan timing (if selected)'
    '                            Fig 3: Slice assignment to volumes'
    '                            Fig 4: cutout interval of logfile for'
    '                            regressor creation (including scan timing'
    '                            and raw phys data)'
    '                            Fig 5: heart beat/breathing statistics & outlier;'
    '                            Fig 6: cardiac phase data of all slices'
    '                            Fig 7: respiratory phase data and'
    '                                   histogram transfer function'
    '                            Fig 8: time course of all sampled RETROICOR'
    '                                   regressors'
    '                            Fig 9: final multiple_regressors matrix'
     
};
verbose.val    = {level fig_output_file use_tabs};




%==========================================================================
%% Structure physio Assemblance
%==========================================================================


%--------------------------------------------------------------------------
% physio
%--------------------------------------------------------------------------
physio      = cfg_exbranch;
physio.tag  = 'physio';
physio.name = 'TAPAS PhysIO Toolbox';
physio.val  = {files sqpar model thresh verbose};
physio.help = {'...'};
physio.prog = @run_physio;
physio.vout = @vout_physio;


%==========================================================================
% function out = run_physio(job)
%==========================================================================
function out = run_physio(job)

job.verbose.fig_handles = [];
[~, R] = tapas_physio_main_create_regressors(job.log_files,  ...
    job.sqpar, job.model, job.thresh, job.verbose);

out.physnoisereg = job.model.output_multiple_regressors;
out.R = R;


%==========================================================================
% function dep = vout_physio(job)
%==========================================================================
function dep = vout_physio(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'physiological noise regressors file';
dep(1).src_output = substruct('.','physnoisereg');
dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
