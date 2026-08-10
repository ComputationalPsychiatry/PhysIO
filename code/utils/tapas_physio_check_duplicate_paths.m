function pathPhysioRoots = tapas_physio_check_duplicate_paths(pathPhysioRoots)
% Checks that only one PhysIO installation is active on the MATLAB path
%
%   pathPhysioRoots = tapas_physio_check_duplicate_paths()
%
% If multiple installations are found, an error lists them and provides a
% ready-to-paste command for keeping each installation in turn.
%
% IN
%   pathPhysioRoots  optional list of installation roots; primarily useful
%                    for testing. By default, all tapas_physio_init files
%                    on the MATLAB path are used.
%
% OUT
%   pathPhysioRoots  canonical paths of the distinct PhysIO installations

% Copyright (C) 2026 TNU, Institute for Biomedical Engineering,
%                    University of Zurich and ETH Zurich.

if nargin < 1
    pathInitFiles = which('tapas_physio_init', '-all');
    if ischar(pathInitFiles)
        pathInitFiles = {pathInitFiles};
    end
    pathPhysioRoots = cellfun(@fileparts, pathInitFiles, ...
        'UniformOutput', false);
elseif ischar(pathPhysioRoots)
    pathPhysioRoots = {pathPhysioRoots};
end

pathPhysioRoots = cellfun(@get_canonical_path, pathPhysioRoots, ...
    'UniformOutput', false);

pathComparison = pathPhysioRoots;
if ispc
    pathComparison = cellfun(@lower, pathComparison, 'UniformOutput', false);
end
[~, indexUnique] = unique(pathComparison, 'stable');
pathPhysioRoots = pathPhysioRoots(indexUnique);

if numel(pathPhysioRoots) <= 1
    return;
end

message = sprintf(['Multiple PhysIO installations are active on the MATLAB path:\n\n' ...
    '%s\nChoose one installation to keep.\n'], format_path_list(pathPhysioRoots));

for iRoot = 1:numel(pathPhysioRoots)
    command = sprintf('cd(''%s''); ', escape_quotes(pathPhysioRoots{iRoot}));
    indexRemove = setdiff(1:numel(pathPhysioRoots), iRoot, 'stable');
    for iRemove = indexRemove
        command = [command, sprintf('rmpath(genpath(''%s'')); ', ...
            escape_quotes(pathPhysioRoots{iRemove}))]; %#ok<AGROW>
    end
    command = [command, 'clear functions; rehash; tapas_physio_init();']; %#ok<AGROW>
    message = [message, sprintf('\nTo keep installation %d, paste:\n%s\n', ...
        iRoot, command)]; %#ok<AGROW>
end

error('tapas:physio:MultipleInstallations', '%s', message);


function pathCanonical = get_canonical_path(pathInput)
% Resolve relative components and filesystem links where possible.
try
    pathFile = java.io.File(pathInput);
    pathCanonical = char(pathFile.getCanonicalPath());
catch
    pathCanonical = char(pathInput);
end


function pathList = format_path_list(pathPhysioRoots)
pathList = '';
for iRoot = 1:numel(pathPhysioRoots)
    pathList = [pathList, sprintf('%d. %s\n', iRoot, ...
        pathPhysioRoots{iRoot})]; %#ok<AGROW>
end


function value = escape_quotes(value)
value = strrep(value, '''', '''''');
