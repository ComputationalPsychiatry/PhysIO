function tests = tapas_physio_conv_test()
% Tests whether current findpeaks function of Matlab's signal processing
% toolbox delivers same results as previous version used in reference data
%
%    tests = tapas_physio_conv_test()
%
% IN
%
% OUT
%
% EXAMPLE
%   tapas_physio_conv_test
%
%   See also

% Author:   Sam Harrison
% Created:  2019-07-17
% Copyright (C) 2019 TNU, Institute for Biomedical Engineering,
%                    University of Zurich and ETH Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under
% the terms of the GNU General Public License (GPL), version 3. You can
% redistribute it and/or modify it under the terms of the GPL (either
% version 3 or, at your option, any later version). For further details,
% see the file COPYING or <http://www.gnu.org/licenses/>.

tests = functiontests(localfunctions);
end

function test_conv_causal(testCase)
%% Tests causal convolution via an impulse response

impulse  = [0 0 0 0 1 0 0 0 0];
filter   = [1 2 3];
solution = [0 0 0 0 1 2 3 0 0];

verifyEqual(testCase, ...
    tapas_physio_conv(impulse, filter, 'causal', 'zero'), ...
    solution);

end

function test_conv_symmetric(testCase)
%% Tests non-causal convolution via an impulse response

impulse  = [0 0 0 0 1 0 0 0 0];
filter   = [1 2 3];
solution = [0 0 0 1 2 3 0 0 0];

verifyEqual(testCase, ...
    tapas_physio_conv(impulse, filter, 'symmetric', 'zero'), ...
    solution);

end