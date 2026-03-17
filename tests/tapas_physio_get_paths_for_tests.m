function [pathExamples, pathTestReferenceResults] = tapas_physio_get_paths_for_tests( ...
    doUseZenodoPaths, doVerifyPath, doDownloadData)
% Returns path to PhysIO test reference results folder, based on
% location of public PhysIO directory, and enables optional download of reference results
% if non-existing
%
%   pathTestReferenceResults = tapas_physio_get_path_test_reference_results(pathPhysioPublic)
%
% IN
%   doUseZenodoPaths    true (1, [default]) or false
%                       if true, use default locations where examples/test
%                       reference results are downloaded to
%                       (PhysIO/examples and PhysIO/test-reference-results)
%                       if false, try to use separate repository locations
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
%   tapas_physio_get_path_test_reference_results(0, 1, 1)
%
%   See also tapas_physio_download_example_data

% Author:   Lars Kasper
% Created:  2026-02-224
% Copyright (C) 2026 TNU, Institute for Biomedical Engineering,
%                    University of Zurich and ETH Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under
% the terms of the GNU General Public License (GPL), version 3. You can
% redistribute it and/or modify it under the terms of the GPL (either
% version 3 or, at your option, any later version). For further details,
% see the file COPYING or <http://www.gnu.org/licenses/>.


if doUseZenodoPaths
    pathPhysioPublic = []; % use default
else
    % leaf directory will be removed
    pathPhysioPublic = 'C:\Users\kasperla\OneDrive - University of Toronto\Documents\Personal\Projects\PhysIO\code\TMP';
end
pathExamples = tapas_physio_get_path_examples(pathPhysioPublic, doVerifyPath, doDownloadData);
pathTestReferenceResults = tapas_physio_get_path_test_reference_results(...
    pathPhysioPublic, doVerifyPath, doDownloadData);
