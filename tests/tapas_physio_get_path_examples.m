function pathExamples = tapas_physio_get_path_examples(pathPhysioPublic, ...
    doVerifyPath, doDownloadData)
% Returns GitLab-internal or TAPAS-public PhysIO Examples folder, based on
% location of public PhysIO directory
%
%   pathExamples = tapas_physio_get_path_examples(pathPhysioPublic)
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
% OUT
%
% EXAMPLE
%   tapas_physio_get_path_examples('tapas/PhysIO')
%
%   See also tapas_physio_download_example_data

% Author:   Lars Kasper
% Created:  2022-09-05
% Copyright (C) 2022 TNU, Institute for Biomedical Engineering,
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

% try PhysIO examples from separate Gitlab repository first (deployment
% mode physio-dev)
% then: try old separate Gitlab repository via submodule in PhysIO Gitlab
% then: try download target from Zenodo download
possibleExamplePaths = {
    fullfile(pathPhysioPublic, '..', 'physio-examples')
    fullfile(pathPhysioPublic, '..', 'examples')
    fullfile(pathPhysioPublic, 'examples')
    };

if doVerifyPath
    iPath = 0;
    haveFoundPath = false;
    while ~haveFoundPath && iPath < numel(possibleExamplePaths)
        iPath = iPath + 1;
        pathExamples = possibleExamplePaths{iPath};
        haveFoundPath = isfolder(fullfile(pathExamples, 'BIDS'));
    end

else % use canonical Zenodo download path for examples
    pathExamples = possibleExamplePaths{end};
end

pathExamples = tapas_physio_simplify_path(pathExamples);

if doDownloadData && (doVerifyPath && ~haveFoundPath)
    pathExamples = tapas_physio_download_example_data();
end

if (doVerifyPath && ~haveFoundPath)
    if doDownloadData
        pathExamples = tapas_physio_download_test_reference_results();
    else
        physio = tapas_physio_new();
        tapas_physio_log('No PhysIO examples data found. Please download via tapas_physio_download_example_data()', physio.verbose, 2);
    end
end

