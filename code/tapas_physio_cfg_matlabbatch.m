function physio = tapas_physio_cfg_matlabbatch
% Lars Kasper, March 2013
%
% Copyright (C) 2013, Institute for Biomedical Engineering, ETH/Uni Zurich.
%
% This file is part of the T, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id$




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
cardiac.tag     = 'log_cardiac';
cardiac.name    = 'log_cardiac';
cardiac.val     = {{'/Users/kasperla/Documents/code/matlab/smoothing_trunk/tSNR_fMRI_SPM/CheckPhysRETROICOR/PhysIOToolbox/examples/Philips/ECG3T/SCANPHYSLOG.log'}};
cardiac.help    = {'...'};
cardiac.filter  = 'any';
cardiac.ufilter = '.*';
cardiac.num     = [1 1];

%--------------------------------------------------------------------------
% respiration (filename)
%--------------------------------------------------------------------------
respiration         = cfg_files;
respiration.tag     = 'log_respiration';
respiration.name    = 'log_respiration';
respiration.val     = {{'/Users/kasperla/Documents/code/matlab/smoothing_trunk/tSNR_fMRI_SPM/CheckPhysRETROICOR/PhysIOToolbox/examples/Philips/ECG3T/SCANPHYSLOG.log'}};
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
sampling_interval.help    = {'sampling interval of phys log files (in seconds)'};
sampling_interval.strtype = 'e';
sampling_interval.num     = [Inf Inf];
sampling_interval.val     = {2e-3};

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
% Nprep
%--------------------------------------------------------------------------
Nprep         = cfg_entry;
Nprep.tag     = 'Nprep';
Nprep.name    = 'Nprep';
Nprep.help    = {'...'};
Nprep.strtype = 'e';
Nprep.num     = [Inf Inf];
Nprep.val     = {[]};

%--------------------------------------------------------------------------
% time_slice_to_slice
%--------------------------------------------------------------------------
time_slice_to_slice         = cfg_entry;
time_slice_to_slice.tag     = 'time_slice_to_slice';
time_slice_to_slice.name    = 'time_slice_to_slice';
time_slice_to_slice.help    = {'...'};
time_slice_to_slice.strtype = 'e';
time_slice_to_slice.num     = [Inf Inf];
time_slice_to_slice.val     = {[]};

%--------------------------------------------------------------------------
% onset_slice
%--------------------------------------------------------------------------
onset_slice         = cfg_entry;
onset_slice.tag     = 'onset_slice';
onset_slice.name    = 'onset_slice';
onset_slice.help    = {'...'};
onset_slice.strtype = 'e';
onset_slice.num     = [Inf Inf];
onset_slice.val     = {19};

%--------------------------------------------------------------------------
% Nscans
%--------------------------------------------------------------------------
Nscans         = cfg_entry;
Nscans.tag     = 'Nscans';
Nscans.name    = 'Nscans';
Nscans.help    = {'...'};
Nscans.strtype = 'e';
Nscans.num     = [Inf Inf];
Nscans.val     = {495};

%--------------------------------------------------------------------------
% Ndummies
%--------------------------------------------------------------------------
Ndummies         = cfg_entry;
Ndummies.tag     = 'Ndummies';
Ndummies.name    = 'Ndummies';
Ndummies.help    = {'...'};
Ndummies.strtype = 'e';
Ndummies.num     = [Inf Inf];
Nscans.val     = {3};

%--------------------------------------------------------------------------
% TR
%--------------------------------------------------------------------------
TR         = cfg_entry;
TR.tag     = 'TR';
TR.name    = 'TR';
TR.help    = {'Repetition time in seconds'};
TR.strtype = 'e';
TR.num     = [Inf Inf];
TR.val     = {2.5};

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
Nscans.val     = {[]};


%--------------------------------------------------------------------------
% Nslices
%--------------------------------------------------------------------------
Nslices         = cfg_entry;
Nslices.tag     = 'Nslices';
Nslices.name    = 'Nslices';
Nslices.help    = {'Number of slices in one volume'};
Nslices.strtype = 'e';
Nslices.num     = [Inf Inf];

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
orthog.tag    = 'orthog';
orthog.name   = 'orthog';
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
input_other_multiple_regressors.filter  = 'mat';
input_other_multiple_regressors.ufilter = '.*';
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


%--------------------------------------------------------------------------
% manual_peak_select
%--------------------------------------------------------------------------
manual_peak_select        = cfg_menu;
manual_peak_select.tag    = 'manual_peak_select';
manual_peak_select.name   = 'manual_peak_select';
manual_peak_select.help   = {'...'};
manual_peak_select.labels = {'true' 'false'};
manual_peak_select.values = {true false};
manual_peak_select.val    = {true};

%--------------------------------------------------------------------------
% kRpeakfile
%--------------------------------------------------------------------------
kRpeakfile         = cfg_files;
kRpeakfile.tag     = 'kRpeakfile';
kRpeakfile.name    = 'kRpeakfile';
kRpeakfile.val     = {{''}};
kRpeakfile.help    = {'...'};
kRpeakfile.filter  = 'any';
kRpeakfile.ufilter = '.*';
kRpeakfile.num     = [0 1];

%--------------------------------------------------------------------------
% min
%--------------------------------------------------------------------------
min       = cfg_entry;
min.tag     = 'min';
min.name    = 'min';
min.help    = {'...'};
min.strtype = 'e';
min.num     = [Inf Inf];

%--------------------------------------------------------------------------
% modality
%--------------------------------------------------------------------------
modality        = cfg_menu;
modality.tag    = 'modality';
modality.name   = 'modality';
modality.help   = {'...'};
modality.labels = {'ECG', 'OXY/PPU'};
modality.values = {'ECG', 'PPU'};
modality.val    = {'ECG'};

%--------------------------------------------------------------------------
% cardiac
%--------------------------------------------------------------------------
cardiac      = cfg_branch;
cardiac.tag  = 'cardiac';
cardiac.name = 'cardiac';
cardiac.val  = {modality min kRpeakfile manual_peak_select};
cardiac.help = {'...'};

%--------------------------------------------------------------------------
% grad_direction
%--------------------------------------------------------------------------
grad_direction        = cfg_menu;
grad_direction.tag    = 'grad_direction';
grad_direction.name   = 'grad_direction';
grad_direction.help   = {'...'};
grad_direction.labels = {'use nominal timing' 'x' 'y' 'z'};
grad_direction.values = {[] 'x' 'y' 'z'};
grad_direction.val    = {[]};

%--------------------------------------------------------------------------
% vol_spacing
%--------------------------------------------------------------------------
vol_spacing         = cfg_entry;
vol_spacing.tag     = 'vol_spacing';
vol_spacing.name    = 'vol_spacing';
vol_spacing.help    = {'...'};
vol_spacing.strtype = 'e';
vol_spacing.num     = [Inf Inf];
vol_spacing.val     = {[]};

%--------------------------------------------------------------------------
% vol
%--------------------------------------------------------------------------
vol         = cfg_entry;
vol.tag     = 'vol';
vol.name    = 'vol';
vol.help    = {'...'};
vol.strtype = 'e';
vol.num     = [Inf Inf];
vol.val     = {[]};

%--------------------------------------------------------------------------
% slice
%--------------------------------------------------------------------------
slice         = cfg_entry;
slice.tag     = 'slice';
slice.name    = 'slice';
slice.help    = {'...'};
slice.strtype = 'e';
slice.num     = [Inf Inf];
slice.val     = {[]};

%--------------------------------------------------------------------------
% zero
%--------------------------------------------------------------------------
zero         = cfg_entry;
zero.tag     = 'zero';
zero.name    = 'zero';
zero.help    = {'...'};
zero.strtype = 'e';
zero.num     = [Inf Inf];
zero.val     = {[]};

%--------------------------------------------------------------------------
% scan_timing
%--------------------------------------------------------------------------
scan_timing      = cfg_branch;
scan_timing.tag  = 'scan_timing';
scan_timing.name = 'scan_timing';
scan_timing.val  = {grad_direction zero slice vol vol_spacing};
scan_timing.help = {'...'};

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
% verbose
%--------------------------------------------------------------------------
level         = cfg_entry;
level.tag     = 'zero';
level.name    = 'zero';
level.help    = {'...'};
level.strtype = 'e';
level.num     = [Inf Inf];
level.val     = {[]};


%--------------------------------------------------------------------------
% verbose
%--------------------------------------------------------------------------
verbose        = cfg_branch;
verbose.tag    = 'verbose';
verbose.name   = 'Display informative figures';
verbose.help   = {'...'};
verbose.val    = {level};




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
[~, R] = tapas_physio_main_create_regressors(job.files,  ...
    job.sqpar, job.model, job.thresh, job.verbose);

out.physnoisereg = job.files.output_multiple_regressors;
out.R = R;


%==========================================================================
% function dep = vout_physio(job)
%==========================================================================
function dep = vout_physio(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'physiological noise regressors file';
dep(1).src_output = substruct('.','physnoisereg');
dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
