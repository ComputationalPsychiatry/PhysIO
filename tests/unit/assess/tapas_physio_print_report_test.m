classdef tapas_physio_print_report_test < matlab.unittest.TestCase
    % Tests extension-aware contrast report printing.

    methods (TestMethodSetup)
        function createGraphicsFigure(testCase)
            graphicsFigure = figure('Visible', 'off', 'Tag', 'Graphics');
            graphicsAxes = axes('Parent', graphicsFigure);
            plot(graphicsAxes, 1:3);
            testCase.addTeardown(@() close(graphicsFigure));
        end
    end

    methods (Test)
        function testExtensionlessReportDefaultsToPdf(testCase)
            outputFolder = testCase.createTemporaryFolder();
            requestedFile = fullfile(outputFolder, 'report');

            fileReport = tapas_physio_print_report(requestedFile);

            testCase.verifyEqual(fileReport, [requestedFile '.pdf']);
            testCase.verifyTrue(isfile(fileReport));
        end

        function testPdfPagesAreAppended(testCase)
            outputFolder = testCase.createTemporaryFolder();
            requestedFile = fullfile(outputFolder, 'report.pdf');
            tapas_physio_print_report(requestedFile);
            initialFile = dir(requestedFile);

            fileReport = tapas_physio_print_report(requestedFile);
            appendedFile = dir(requestedFile);

            testCase.verifyEqual(fileReport, requestedFile);
            testCase.verifyGreaterThan(appendedFile.bytes, initialFile.bytes);
        end

        function testUnsupportedExtensionErrors(testCase)
            outputFolder = testCase.createTemporaryFolder();
            requestedFile = fullfile(outputFolder, 'report.unsupported');

            testCase.verifyError( ...
                @() tapas_physio_print_report(requestedFile), ...
                'tapas_physio:UnsupportedPrintFormat');
        end
    end
end
