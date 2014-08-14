function [correlation,x,y] = tapas_physio_corrcoef12(x,y, isZtransformed)
% computes correlation coefficient (i.e. entry (1,2) of correlation matrix) 
% quickly between two time series
% 
%
%   [correlation,x,y] = tapas_physio_corrcoef12(x,y, isZtransformed)
%
%
% IN
%   x               [nSamples,1] column vector of samples
%   y               [nSamples,1] column vector of samples
%   isZtransformed [1,2] vector stating whether x,y or both are
%                   z-transformed, i.e. mean-corrected and divided by their
%                   standard deviation
%                   example: 
%                   isZtransformed = [1,0] means that x is z-transformed,
%                   by y is not, i.e. (y-mean(y))/std(y) will be computed
% OUT
%   correlation     correlation(1,2) of correlation = corrcoef(x,y)
%   x               z-transformed x
%   y               z-transformed y
% EXAMPLE
%   tapas_physio_corrcoef12
%
%   See also
%
% Author: Lars Kasper
% Created: 2014-08-14
% Copyright (C) 2014 TNU, Institute for Biomedical Engineering, University of Zurich and ETH Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id$
doUseSlow = false;

if nargin < 3
    isZtransformed = [0 0];
end


if doUseSlow
    correlation = corrcoef(x,y);
    correlation = correlation(1,2);
else %fast, using shortcuts...
    %C(i,j)/SQRT(C(i,i)*C(j,j)).
    
    x = x(:);
    y = y(:);
    nSamples = numel(x);

    if ~isZtransformed(1)
       x = x - mean(x); x = x./sqrt(x'*x/(nSamples-1));
    end
    if ~isZtransformed(2)
        y = y - mean(y); y = y./sqrt(y'*y/(nSamples-1));
    end
    correlation = x'*y;
end