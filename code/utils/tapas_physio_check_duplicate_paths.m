function pathPhysioRoots = tapas_physio_check_duplicate_paths(pathPhysioRoots)
% Checks that only one PhysIO installation is active on the MATLAB path
%
%   pathPhysioRoots = tapas_physio_check_duplicate_paths()
%
% If multiple installations are found, recovery instructions are printed
% before an error is raised. Each instruction line is a complete MATLAB
% command, so the whole block can be pasted into the Command Window.
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

instructions = sprintf(['Multiple PhysIO installations are active on the MATLAB path:\n\n' ...
    '%s\nChoose one installation to keep.\n'], format_path_list(pathPhysioRoots));

for iRoot = 1:numel(pathPhysioRoots)
    command = sprintf('cd(''%s'')\n', escape_quotes(pathPhysioRoots{iRoot}));
    indexRemove = setdiff(1:numel(pathPhysioRoots), iRoot, 'stable');
    for iRemove = indexRemove
        command = [command, sprintf('rmpath(genpath(''%s''))\n', ...
            escape_quotes(pathPhysioRoots{iRemove}))]; %#ok<AGROW>
    end
    command = [command, sprintf('clear functions\nrehash\ntapas_physio_init()\n')]; %#ok<AGROW>
    instructions = [instructions, sprintf('\nTo keep installation %d, paste:\n%s', ...
        iRoot, command)]; %#ok<AGROW>
end

% ERROR reformats long messages and may insert line breaks inside quoted
% paths. FPRINTF preserves the commands exactly as generated.
fprintf(2, '\n%s\n', instructions);
% Supply verbose explicitly: the logger's default calls tapas_physio_new,
% which checks duplicate installations again and would recurse here.
verbose = struct('level', 0, 'process_log', {{instructions}});
tapas_physio_log( ...
    'Multiple PhysIO installations are active. See recovery commands above.', ...
    verbose, 2);


function pathCanonical = get_canonical_path(pathInput)
% Resolve relative components and filesystem links where possible.
pathCanonical = char(pathInput);
if exist('resolveFilePath', 'builtin') || exist('resolveFilePath', 'file')
    pathResolved = resolveFilePath(pathCanonical, 'EntryType', 'folder', ...
        'ResolveSymbolicLinks', true);
    if ~ismissing(pathResolved)
        pathCanonical = char(pathResolved);
        return;
    end
end

% Older releases can canonicalize with Java, but -nojvm sessions cannot.
% Avoid even attempting Java construction without a JVM: repeated caught
% JavaNotFound errors can break subsequent path updates in R2026b.
if usejava('jvm')
    try
        pathFile = java.io.File(pathCanonical);
        pathCanonical = char(pathFile.getCanonicalPath());
        return;
    catch
        % Fall back to MATLAB's absolute folder name below.
    end
end
% FILEATTRIB returns an absolute name for an existing folder, including an
% empty one; unlike the resolvers above, it does not guarantee link resolution.
[success, attributes] = fileattrib(pathCanonical);
if success && isscalar(attributes) && attributes.directory
    pathCanonical = attributes.Name;
end


function pathList = format_path_list(pathPhysioRoots)
pathList = '';
for iRoot = 1:numel(pathPhysioRoots)
    pathList = [pathList, sprintf('%d. %s\n', iRoot, ...
        pathPhysioRoots{iRoot})]; %#ok<AGROW>
end


function value = escape_quotes(value)
value = strrep(value, '''', '''''');
