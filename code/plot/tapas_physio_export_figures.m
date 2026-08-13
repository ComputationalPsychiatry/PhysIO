function fileName = tapas_physio_export_figures( ...
        figHandles, fileName, appendToExisting)
% Export figure handles to a multipage PDF or legacy PostScript file.
%
%   fileName = tapas_physio_export_figures( ...
%       figHandles, fileName, appendToExisting)
%
% PDF export uses exportgraphics on MATLAB R2020a and newer. On those
% releases, a requested PostScript file is converted to PDF. Older MATLAB
% releases use PostScript because multipage PDF export is unavailable.

if nargin < 3
    appendToExisting = false;
end

[pathName, baseName, extension] = fileparts(fileName);
switch lower(extension)
    case '.ps'
        if ~verLessThan('matlab', '9.8')
            warning('tapas:physio:PostScriptNotSupported', ...
                ['PostScript export is no longer supported by this ' ...
                'MATLAB release. Saving the report as PDF instead.']);
            fileName = fullfile(pathName, [baseName '.pdf']);
            export_to_pdf(figHandles, fileName, appendToExisting);
        else
            export_to_postscript(figHandles, fileName, appendToExisting);
        end
    case '.pdf'
        if ~verLessThan('matlab', '9.8')
            export_to_pdf(figHandles, fileName, appendToExisting);
        else
            warning('tapas:physio:PdfNotSupported', ...
                ['Multipage PDF export requires MATLAB R2020a or ' ...
                'newer. Saving the report as PostScript instead.']);
            fileName = fullfile(pathName, [baseName '.ps']);
            export_to_postscript(figHandles, fileName, appendToExisting);
        end
    otherwise
        error('tapas_physio:UnsupportedMultipageFormat', ...
            'Multipage figure export does not support extension: %s', ...
            extension);
end
end

function export_to_pdf(figHandles, fileName, appendToExisting)
for k = 1:numel(figHandles)
    doAppend = appendToExisting || k > 1;
    exportgraphics(figHandles(k), fileName, 'Append', doAppend);
end
end

function export_to_postscript(figHandles, fileName, appendToExisting)
if ~appendToExisting && isfile(fileName)
    delete(fileName);
end
try % Level 2 PostScript
    for k = 1:numel(figHandles)
        print(figHandles(k), '-dpsc2', '-append', fileName); %#ok<PRINTPS4>
    end
catch
    if isfile(fileName)
        delete(fileName);
    end
    for k = 1:numel(figHandles)
        print(figHandles(k), '-dpsc', '-append', fileName); %#ok<PRINTPS2>
    end
end
end
