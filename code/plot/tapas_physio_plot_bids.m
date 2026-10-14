function [verbose, ons_secs, sqpar] = tapas_physio_plot_bids(filename, sqpar, verbose, spulse)
% Plots BIDS recordings using the existing PhysIO review figures
%
%   [verbose, ons_secs, sqpar] = tapas_physio_plot_bids(filename, ...
%       sqpar, verbose, spulse)
%
% IN
%   filename    .json, .tsv or .tsv.gz file; the matching data file and
%               JSON sidecar are inferred from the filename
%   sqpar       optional scan parameters (default: empty structure)
%   .Ndummies   number of dummy volumes (default: 0)
%   .Nslices    number of slices per volume (default: 1)
%   verbose     optional figure settings from tapas_physio_new
%   spulse      optional measured slice-event times [NslicesTotal, 1],
%               in seconds relative to the first volume; required when
%               sqpar.Nslices > 1. For one slice, defaults to volume onsets
%
% OUT
%   verbose     figure settings with the two figures added to .fig_handles
%   ons_secs    read-in physiological traces and events passed to the plots
%   .raw        identical copy of the read-in data for the raw-data overlay
%   .fr         same as .r; no respiratory filtering is performed
%   sqpar       scan parameters with defaults filled in
%               Recorded event counts are available as
%               numel(ons_secs.svolpulse) and numel(ons_secs.cpulse).
%
% EXAMPLE
%   [verbose, ons_secs] = tapas_physio_plot_bids('sub-01_physio.json');
%
% No filtering, cropping or additional heartbeat detection is performed.
% Empty optional arguments select defaults. Slice timing cannot be recovered
% from volume triggers alone. Cardiac diagnostics use the default posthoc
% pulse-selection thresholds from tapas_physio_new.
%
%   See also tapas_physio_read_physlogfiles_bids
%       tapas_physio_plot_cropped_phys_to_acqwindow
%       tapas_physio_plot_raw_physdata_diagnostics

% Copyright (C) 2026, University of Zurich and ETH Zurich.
% This file is part of PhysIO, released under the GNU GPL version 3 or later.

defaults = tapas_physio_new();
if nargin < 2 || isempty(sqpar), sqpar = struct(); end
if nargin < 3 || isempty(verbose), verbose = defaults.verbose; end
if nargin < 4, spulse = []; end
% The BIDS physiological JSON/TSV files do not specify the number of dummy
% volumes or slices per volume, or the timing of individual slice events.
% The acquisition-window plot requires these fields: unless supplied by the
% caller, assume no dummy volumes and one slice per volume. Below, volume
% onsets also serve as slice events for this one-slice plotting fallback;
% they do not represent measured slice timing. For multiple slices, the
% caller must supply both sqpar.Nslices and the measured spulse times.
if ~isfield(sqpar, 'Ndummies') || isempty(sqpar.Ndummies)
    sqpar.Ndummies = 0;
end
if ~isfield(sqpar, 'Nslices') || isempty(sqpar.Nslices)
    sqpar.Nslices = 1;
end
filename = char(filename);
base = regexprep(filename, '\.(json|tsv(\.gz)?)$', '', 'ignorecase');
if strcmp(base, filename)
    verbose = tapas_physio_log('Expected .json, .tsv or .tsv.gz.', verbose, 2);
end
fileJson = [base '.json'];
if endsWith(lower(filename), '.json')
    filename = [base '.tsv.gz'];
    if ~isfile(filename), filename = [base '.tsv']; end
end
if ~isfile(filename) || ~isfile(fileJson)
    verbose = tapas_physio_log(...
        'BIDS data and JSON sidecar are required.', verbose, 2);
end
log_files = defaults.log_files;
log_files.cardiac = filename;
log_files.respiration = filename;
log_files.scan_timing = fileJson;
log_files.relative_start_acquisition = [];
readerVerbose = verbose;
readerVerbose.level = 0;
[ons_secs.c, ons_secs.r, ons_secs.t, ons_secs.cpulse, ...
    ons_secs.acq_codes] = tapas_physio_read_physlogfiles_bids(...
    log_files, 'ECG', readerVerbose);
ons_secs.svolpulse = ons_secs.t(bitand(uint32(ons_secs.acq_codes), 8) ~= 0);
ons_secs.svolpulse = ons_secs.svolpulse(:);
if isempty(spulse)
    if sqpar.Nslices ~= 1
        verbose = tapas_physio_log(...
            'Supply measured spulse when sqpar.Nslices is greater than one.', ...
            verbose, 2);
    end
    spulse = ons_secs.svolpulse;
end
ons_secs.spulse = spulse(:);
if sqpar.Ndummies < 0 || fix(sqpar.Ndummies) ~= sqpar.Ndummies || ...
        sqpar.Ndummies > numel(ons_secs.svolpulse) || ...
        sqpar.Ndummies * sqpar.Nslices > numel(ons_secs.spulse)
    verbose = tapas_physio_log(...
        'Dummy scans exceed available events.', verbose, 2);
end
% Report recorded events before plotting; volume markers do not necessarily
% establish the number of acquired volumes if triggers are missing.
verbose = tapas_physio_log(sprintf(...
    'BIDS recording: %d samples, %d volume events, %d cardiac pulse events', ...
    numel(ons_secs.t), numel(ons_secs.svolpulse), numel(ons_secs.cpulse)), ...
    verbose, 0);
ons_secs.fr = ons_secs.r;
ons_secs.t_start = ons_secs.t(1);
ons_secs.raw = ons_secs;
fh = tapas_physio_plot_cropped_phys_to_acqwindow(ons_secs, sqpar, verbose);
set(fh, 'Name', 'BIDS: Recorded physiological signals and events');
verbose.fig_handles(end+1) = fh;
% The legacy diagnostics plot requires a waveform even for pulse-only data.
cDiagnostics = ons_secs.c;
if isempty(cDiagnostics)
    cDiagnostics = double(ismember(ons_secs.t, ons_secs.cpulse));
end
verbose.level = max(1, verbose.level);
verbose = tapas_physio_plot_raw_physdata_diagnostics(ons_secs.cpulse, ...
    ons_secs.r, defaults.preproc.cardiac.posthoc_cpulse_select, ...
    verbose, ons_secs.t, cDiagnostics);
end
