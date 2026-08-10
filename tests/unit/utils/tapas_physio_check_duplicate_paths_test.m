function tests = tapas_physio_check_duplicate_paths_test()
% Unit tests for tapas_physio_check_duplicate_paths

tests = functiontests(localfunctions);


function test_single_installation_is_accepted(testCase)
pathRoot = tempname;
mkdir(pathRoot);
cleanup = onCleanup(@() rmdir(pathRoot));

pathActual = tapas_physio_check_duplicate_paths({pathRoot});

verifyEqual(testCase, pathActual, {get_canonical_path(pathRoot)});


function test_duplicate_installations_report_recovery_commands(testCase)
pathRoot1 = tempname;
pathRoot2 = [tempname, '''s copy'];
mkdir(pathRoot1);
mkdir(pathRoot2);
cleanup = onCleanup(@() remove_test_folders(pathRoot1, pathRoot2));

exception = capture_exception(@() tapas_physio_check_duplicate_paths( ...
    {pathRoot1, pathRoot2}));
pathRoot1 = get_canonical_path(pathRoot1);
pathRoot2 = get_canonical_path(pathRoot2);

verifyEqual(testCase, exception.identifier, ...
    'tapas:physio:MultipleInstallations');
verifySubstring(testCase, exception.message, ['1. ', pathRoot1]);
verifySubstring(testCase, exception.message, ['2. ', pathRoot2]);
verifySubstring(testCase, exception.message, sprintf( ...
    'cd(''%s''); rmpath(genpath(''%s'')); clear functions; rehash; tapas_physio_init();', ...
    pathRoot1, strrep(pathRoot2, '''', '''''')));
verifySubstring(testCase, exception.message, sprintf( ...
    'cd(''%s''); rmpath(genpath(''%s'')); clear functions; rehash; tapas_physio_init();', ...
    strrep(pathRoot2, '''', ''''''), pathRoot1));


function exception = capture_exception(functionHandle)
try
    functionHandle();
catch exception
    return;
end
error('Expected tapas_physio_check_duplicate_paths to throw an error.');


function pathCanonical = get_canonical_path(pathInput)
pathFile = java.io.File(pathInput);
pathCanonical = char(pathFile.getCanonicalPath());


function remove_test_folders(varargin)
for iFolder = 1:nargin
    if exist(varargin{iFolder}, 'dir')
        rmdir(varargin{iFolder});
    end
end
