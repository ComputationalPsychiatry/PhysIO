function fileReport = tapas_physio_print_report(fileReport)
% Print the SPM Graphics window to an extension-aware report file.
%
% If fileReport has no extension, PDF is used. PDF pages are appended with
% exportgraphics on MATLAB R2020a and newer. Other extensions supported by
% spm_print, including explicit PostScript requests, use SPM's exporter.

[pathReport, nameReport, extensionReport] = fileparts(fileReport);
if isempty(extensionReport)
    extensionReport = '.pdf';
end
fileReport = fullfile(pathReport, [nameReport extensionReport]);

if strcmpi(extensionReport, '.pdf') && verLessThan('matlab', '9.8') == 0
    graphicsFigure = spm_figure('FindWin', 'Graphics');
    if isempty(graphicsFigure)
        error('tapas_physio:GraphicsFigureNotFound', ...
            'The SPM Graphics figure could not be found.');
    end
    fileReport = tapas_physio_export_figures( ...
        graphicsFigure, fileReport, isfile(fileReport));
    return;
end

printFormat = extensionReport(2:end);
printOptions = spm_print('format', printFormat);
if isempty(printOptions)
    error('tapas_physio:UnsupportedPrintFormat', ...
        'Unsupported report file extension: %s', extensionReport);
end
spm_print(fileReport, 'Graphics', printOptions);
end
