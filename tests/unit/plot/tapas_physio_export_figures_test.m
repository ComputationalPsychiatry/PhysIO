classdef tapas_physio_export_figures_test < matlab.unittest.TestCase
    % Tests the shared multipage figure exporter.

    methods (Test)
        function testReplacesExistingPdfByDefault(testCase)
            outputFolder = testCase.createTemporaryFolder();
            fileName = fullfile(outputFolder, 'report.pdf');
            firstFigure = create_figure(1:3);
            secondFigure = create_figure(3:-1:1);
            testCase.addTeardown(@() close([firstFigure secondFigure]));
            tapas_physio_export_figures(firstFigure, fileName);
            firstFile = dir(fileName);

            tapas_physio_export_figures(secondFigure, fileName);
            replacedFile = dir(fileName);

            testCase.verifyTrue(isfile(fileName));
            testCase.verifyLessThanOrEqual( ...
                replacedFile.bytes, 2 * firstFile.bytes);
        end

        function testAppendsToExistingPdfWhenRequested(testCase)
            outputFolder = testCase.createTemporaryFolder();
            fileName = fullfile(outputFolder, 'report.pdf');
            graphicsFigure = create_figure(1:3);
            testCase.addTeardown(@() close(graphicsFigure));
            tapas_physio_export_figures(graphicsFigure, fileName);
            initialFile = dir(fileName);

            tapas_physio_export_figures(graphicsFigure, fileName, true);
            appendedFile = dir(fileName);

            testCase.verifyGreaterThan(appendedFile.bytes, initialFile.bytes);
        end
    end
end

function graphicsFigure = create_figure(plotData)
graphicsFigure = figure('Visible', 'off');
graphicsAxes = axes('Parent', graphicsFigure);
plot(graphicsAxes, plotData);
end
