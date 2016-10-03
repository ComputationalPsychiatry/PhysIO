function [tSnrImage, fileTsnr] = tapas_physio_compute_tsnr_spm(SPM, iC, doSave)
% Computes temporal SNR image after correcting for a contrast 
% from SPM general linear model
%
%   tSnrImage = compute_tsnr_spm(input)
%
% IN
%   SPM     SPM variable (in SPM.mat, or file name) after parameter and 
%           contrast estimation
%   iC      Contrast of interest for which tSNR is computed
%
%           iC = 0 computes tSNR of the raw image time series after
%               pre-whitening and filtering, but without any model
%               confound regressors removed. (default)
%           iC = NaN computes tSNR after removing the full model (including
%               the mean), and is therefore equivalent to sqrt(1/ResMS)
%   doInvert true (default for iC > 0 and ~nan) or false
%           typically, one is interested in the tSNR *after* correcting for
%           the contrast in question. To compute this, one has to look at
%           the residuals of the inverse F-contrast that excludes all but
%           the contrast of interest, i.e. eye(nRegressors) - xcon
%           Thus, this function computes this contrast per default and goes
%           from there determining residuals etc.
%
%   doSave  true (default) or false
%           if true, a file tSNR_con<iC>.nii is created in the same folder as
%           the SPM.mat
%   
% OUT
%   tSnrImage   MrImage holding tSNR when only including regressors in
%               contrast iC in design matrix;
%               i.e. gives the effect of
%                       mean(Xc*bc + e)/std(Xc*bc + e)
%               if Xc hold only regressors of iC
%               
%               NOTE: If you want to estimate the effect of a noise
%               correction using some regressors, the Contrast with index
%               iC should actually contain an F-contrast *EXCLUDING* these
%               regressors, because everything else constitutes the
%               regressors of interest
%   fileTsnr    tsnr_con<iC>.nii
%               path and file name of tSNR image, if doSave was true
%
% EXAMPLE
%   compute_tsnr_spm
%
%   See also spm_write_residuals spm_FcUtil
%
% Author: Lars Kasper
% Created: 2014-12-10
% Copyright (C) 2014 Institute for Biomedical Engineering, ETH/Uni Zurich.
% $Id$

if nargin < 2
    iC = 0;
end

if nargin < 3
    doSave = true;
end

if nargin < 4
    doInvert = true;
end

% load SPM-variable, if filename given
if ~isstruct(SPM)
    % load SPM variable from file
    if iscell(SPM)
        fileSpm = SPM{1};
    else
        fileSpm = SPM;
    end
    load(fileSpm, 'SPM');
    
    % temporary changes to SPM structure saved in sub-dir
    oldDirSpm = SPM.swd;
    newDirSpm = fullfile(SPM.swd, 'tmp'); 
    mkdir(newDirSpm);
    copyfile(fullfile(SPM.swd, '*.nii'), newDirSpm);
    copyfile(fullfile(SPM.swd, 'SPM.mat'), newDirSpm);
    SPM.swd = newDirSpm;
end

 iCIn = iC;

 isInvertableContrast = iC > 0 && ~isnan(iC);
 
 if isInvertableContrast && doInvert
     if ~isequal(SPM.xCon(iC).STAT, 'F')
         error('Can only invert F-contrasts');
     end
     
     idxColumnsContrast = find(sum(SPM.xCon(iC).c));
     
     Fc = spm_FcUtil('Set', ['All but: ' SPM.xCon(iC).name], 'F', ...
         'iX0', idxColumnsContrast, SPM.xX.xKXs);
     SPM.xCon(end+1) = Fc;
     SPM = spm_contrasts(SPM,numel(SPM.xCon));
     
     % use this for computation
     iC = numel(SPM.xCon);
 end

%% Write residuals Y - Y0 = Yc + e;
VRes        = spm_write_residuals(SPM, iC);
nVolumes    = numel(VRes);
for iVol = 1:nVolumes
    VRes(iVol).fname = fullfile(SPM.swd, VRes(iVol).fname);
end

fileTsnr = fullfile(oldDirSpm, sprintf('tSNR_con%04d.nii', iCIn));

useMrImage = false;

if useMrImage % use toolbox functionality
    ResImage    = MrImage(VRes(1).fname);
    
    % Create 4D image of contrast-specific "residuals", i.e. Xc*bc + e
    for iVol = 1:nVolumes
        ResImage.append(VRes(iVol).fname);
    end
    
    % compute tSNR = mean(Xc*bc + e)/std(Xc*bc + e)
    tSnrImage = ResImage.mean./ResImage.std;
    if doSave
        tSnrImage.save(fileTsnr);
    end
else
    ResImage = spm_read_vols(VRes);
    meanImage = mean(ResImage, 4);
    stdImage = std(ResImage, 0, 4);
    tSnrImage = meanImage./stdImage;
    VTsnr = VRes(1);
    VTsnr.fname = fileTsnr;
    spm_write_vol(VTsnr, tSnrImage);
end


%% clean up all created residual files and temporary SPM folder
delete(fullfile(newDirSpm, '*'));
rmdir(newDirSpm);