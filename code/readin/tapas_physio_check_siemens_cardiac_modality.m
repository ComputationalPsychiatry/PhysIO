function verbose = tapas_physio_check_siemens_cardiac_modality(fileName, cardiacModality, verbose)
%tapas_physio_check_siemens_cardiac_modality - Check cardiac modality
%   tapas_physio_check_siemens_cardiac_modality(fileName, cardiacModality)
%   raises an error if the Siemens filename indicates PULS but the modality
%   is ECG, or vice versa. Empty or unrecognized filenames are left alone.
%
%   VERBOSE = tapas_physio_check_siemens_cardiac_modality(..., VERBOSE)
%   uses and returns the supplied PhysIO logging settings.
%
%   See also tapas_physio_log

% Copyright (C) 2026 TNU, Institute for Biomedical Engineering,
%                    University of Zurich and ETH Zurich.
% This file is part of PhysIO, released under the GNU GPL, version 3.

if nargin < 3 || isempty(verbose)
    verbose = struct('level', 0, 'process_log', {{}});
end
if isempty(fileName)
    return;
end

[~, name, extension] = fileparts(fileName);
isPulsFile = strcmpi(extension, '.puls') || ...
    ~isempty(regexpi(name, '(^|_)PULS$', 'once'));
isEcgFile = strcmpi(extension, '.ecg') || ...
    ~isempty(regexpi(name, '(^|_)ECG[1-4]?$', 'once'));
isEcg = ~isempty(regexpi(cardiacModality, '^ecg($|_)', 'once'));
isPulseOximetry = any(strcmpi(cardiacModality, {'PULS', 'PPU', 'OXY', 'PMU'}));

if (isPulsFile && isEcg) || (isEcgFile && isPulseOximetry)
    verbose = tapas_physio_log(sprintf(['Cardiac modality ''%s'' conflicts with Siemens file ''%s''. ' ...
        'Set preproc.cardiac.modality to ''PPU'' (or ''PULS''/''OXY'') for PULS files, ' ...
        'or ''ECG'' for ECG files.'], cardiacModality, fileName), verbose, 2);
end
end
