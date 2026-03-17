function zenodoRecord = tapas_physio_get_zenodo_record_for_version(conceptRecid, semVersion)
% Returns Zenodo record fields for specified concept RecId (parent DOI of all
% versions) and PhysIO version
%
%   zenodoRecord = tapas_physio_get_zenodo_record_for_version(conceptRecid, semVersion)
%
% IN
%   conceptRecid        numeric Zenodo concept record identifier that groups
%                       all version-specific PhysIO example-data records
%                       e.g. 15579086 for concept DOI
%                       10.5281/zenodo.15579086
%   semVersion          semantic PhysIO version string to match against the
%                       Zenodo record metadata, e.g. 'v9.2.0'
%                       if empty or omitted, the current PhysIO version from
%                       tapas_physio_version() is used
% OUT
%   zenodoRecord        Zenodo API record structure for the requested
%                       version, including metadata and download links
%                       for example, zenodoRecord.doi_url gives you the full DOI 
%                       for the requested version, e.g., "https://doi.org/10.5281/zenodo.19075176"
%
% EXAMPLE
%   zenodoRecord = tapas_physio_get_zenodo_record_for_version(15579086, 'v9.1.0')
%
%   See also tapas_physio_download_example_data tapas_physio_version
 
% Author:   Johanna Bayer & Lars Kasper
% Created:  2026-03-17
% Copyright (C) 2026 TNU, Institute for Biomedical Engineering,
%                    University of Zurich and ETH Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under
% the terms of the GNU General Public License (GPL), version 3. You can
% redistribute it and/or modify it under the terms of the GPL (either
% version 3 or, at your option, any later version). For further details,
% see the file COPYING or <http://www.gnu.org/licenses/>.

if nargin < 2
    currentRelease = tapas_physio_version();
    semVersion = regexprep(currentRelease, '.*v', 'v');
end

s = webread(sprintf( ...
    'https://zenodo.org/api/records?all_versions=true&q=conceptrecid:%d', ...
    conceptRecid));

hits = s.hits.hits;
if ~iscell(hits), hits = num2cell(hits); end

idx = find(cellfun(@(x) strcmp(x.metadata.version, semVersion), hits), 1);

if isempty(idx)
    error('No Zenodo record found for conceptrecid %d and version %s.', ...
        conceptrecid, version);
end

zenodoRecord = hits{idx};
