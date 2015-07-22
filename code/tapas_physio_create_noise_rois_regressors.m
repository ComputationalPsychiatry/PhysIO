function R_noise_rois = tapas_physio_create_noise_rois_regressors(noise_rois)
% compute physiological regressors by extracting principal components of
% the preprocessed fMRI time series from anatomically defined noise ROIS (e.g. CSF)
%
%   R_noise_rois = tapas_physio_create_noise_rois_regressors(noise_rois)
% 
% Approach similar to the one described as aCompCor:
% Behzadi, Y., Restom, K., Liau, J., Liu, T.T., 2007. 
% A component based noise correction method (CompCor) for BOLD and 
% perfusion based fMRI. NeuroImage 37, 90?101. 
% doi:10.1016/j.neuroimage.2007.04.042
%
% IN
%   physio.model.noise_rois
% OUT
%
% EXAMPLE
%   tapas_physio_create_noise_rois_regressors
%
%   See also
%
% Author: Lars Kasper
% Created: 2015-07-22
% Copyright (C) 2015 TNU, Institute for Biomedical Engineering,
%                    University of Zurich and ETH Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under the terms of the GNU General Public
% License (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id$
