function pathTestReferenceResults = tapas_physio_get_path_test_reference_results(pathPhysioPublic, ...
    doVerifyPath, doDownloadData)
% Returns path to PhysIO test reference results folder, based on
% location of public PhysIO directory, and enables optional download of reference results
% if non-existing
%
%   pathTestReferenceResults = tapas_physio_get_path_test_reference_results(pathPhysioPublic)
%
% IN
%   pathPhysioPublic    location of public PhysIO folder, e.g.,
%                       'tapas/PhysIO'
%                       leave empty if location of PhysIO code should be
%                       determined automatically
%   doVerifyPath        true (1, [default]) or false
%                       if true, a set of generic example paths are tested
%                       and the first one that exists is returned
%                       if false, the standard download path for examples,
%                       as used by tapas_physio_download_example_data(), is
%                       returned
%   doDownloadData      true or false [default]; if true and no valid path 
%                       exists, initiate download of reference data from
%                       web
%
% OUT
%
% EXAMPLE
%   tapas_physio_get_path_test_reference_results()
%
%   See also tapas_physio_download_example_data

% Author:   Lars Kasper
% Created:  2025-07-29
% Copyright (C) 2025 TNU, Institute for Biomedical Engineering,
%                    University of Zurich and ETH Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under
% the terms of the GNU General Public License (GPL), version 3. You can
% redistribute it and/or modify it under the terms of the GPL (either
% version 3 or, at your option, any later version). For further details,
% see the file COPYING or <http://www.gnu.org/licenses/>.

if nargin < 1 || isempty(pathPhysioPublic)
    pathPhysioPublic = tapas_physio_simplify_path(fullfile(fileparts(mfilename('fullpath')), '..'));
end

if nargin < 2
    doVerifyPath = true;
end

if nargin < 3
    doDownloadData = false;
end

% try canonical download locations for test data from current and archived PhysIO release
% conventions
possiblePaths = {
    fullfile(pathPhysioPublic, '..', 'physio-examples', 'TestReferenceResults')
    fullfile(pathPhysioPublic, '..', 'physio-test-reference-results')
    fullfile(pathPhysioPublic, 'test-reference-results')
    };

if doVerifyPath
    iPath = 0;
    haveFoundPath = false;
    while ~haveFoundPath && iPath < numel(possiblePaths)
        iPath = iPath + 1;
        pathTestReferenceResults = possiblePaths{iPath};
        haveFoundPath = isfolder(fullfile(pathTestReferenceResults, 'BIDS'));
    end

    % If no examples folder found, suggest to download them via tapas-function
    if ~isfolder(fullfile(pathTestReferenceResults, 'BIDS'))
        physio = tapas_physio_new();
        tapas_physio_log('No PhysIO examples data found. Please download via tapas_physio_download_example_data()', physio.verbose, 2);
    end

else % use canonical Zenodo download path for examples
    pathTestReferenceResults = possiblePaths{end};
end

pathTestReferenceResults = tapas_physio_simplify_path(pathTestReferenceResults);

if doDownloadData
 
end
