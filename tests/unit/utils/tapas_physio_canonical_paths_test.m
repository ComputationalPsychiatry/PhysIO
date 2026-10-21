classdef tapas_physio_canonical_paths_test < matlab.unittest.TestCase
    methods (Test)
        function relativeAliasesAreOneInstallation(testCase)
            pathRoot = testCase.createTemporaryFolder();
            pathActual = tapas_physio_check_duplicate_paths( ...
                {pathRoot, fullfile(pathRoot, '.')});
            testCase.verifyNumElements(pathActual, 1);
            testCase.verifyTrue(isfolder(pathActual{1}));
        end

        function missingFolderIsPreserved(testCase)
            pathRoot = fullfile(testCase.createTemporaryFolder(), 'missing');
            pathActual = tapas_physio_check_duplicate_paths({pathRoot});
            testCase.verifyEqual(pathActual, {pathRoot});
        end

        function repeatedChecksPreservePathUpdates(testCase)
            pathRoot = testCase.createTemporaryFolder();
            testCase.verifyWarningFree(@() repeat_checks(pathRoot));
        end
    end
end

function repeat_checks(pathRoot)
% The unguarded Java call caused native path failures after ten iterations
% in R2026b without a JVM, even though its MATLAB exception was caught.
for iteration = 1:30
    tapas_physio_check_duplicate_paths({pathRoot});
    path(path);
    drawnow;
end
end
