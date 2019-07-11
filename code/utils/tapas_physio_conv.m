function [w] = tapas_physio_conv(u, v, padding)
% Wrapper around `conv()` for causal convolution
%   Deals with the padding and time offsets
%
%    [w] = tapas_physio_filter_cardiac(u, v)
%
% IN
%   u           Data time series [1-D]
%   v           Convolutional filter time series [1-D]
%               Note this is defined for `t >= 0`, with the first element
%               corresponding to `t=0`.
%   padding     ['mean', 'zero']. Whether to pad with the mean of `u`, or
%               with zeros.
%
% OUT
%   w           Convolved time series [size(u)]

% Author:   Sam Harrison
% Created:  2019-07-11
% Copyright (C) 2019 TNU, Institute for Biomedical Engineering,
%                    University of Zurich and ETH Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under
% the terms of the GNU General Public License (GPL), version 3. You can
% redistribute it and/or modify it under the terms of the GPL (either
% version 3 or, at your option, any later version). For further details,
% see the file COPYING or <http://www.gnu.org/licenses/>.

if ~isvector(u) || ~isvector(v)
    error('tapas_physio_conv: Both inputs must be vectors')
end
if nargin < 3
    padding = 'mean';
end

% Padding: use data mean to reduce transients
switch lower(padding)
    case 'mean'
        pad_val = mean(u);
    case 'zero'
        pad_val = 0.0;
    otherwise
        error('Unrecognised padding argument (%s)', padding)
end
u = [pad_val * ones(length(v)-1, 1); u(:)];

% Apply convolution and select portion of interest
w = conv(u, v, 'valid');

end