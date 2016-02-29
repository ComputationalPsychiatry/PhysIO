function [c, r, t, cpulse, verbose] = tapas_physio_read_physlogfiles_siemens(log_files, ...
    verbose, varargin)
% reads out physiological time series and timing vector for Siemens
% logfiles of peripheral cardiac monitoring (ECG/Breathing Belt or
% pulse oximetry)
%
%   [cpulse, rpulse, t, c] = tapas_physio_read_physlogfiles_siemens(logfile, vendor, cardiac_modality)
%
% IN    log_files
%       .log_cardiac        contains ECG or pulse oximeter time course
%                           for GE: ECGData...
%       .log_respiration    contains breathing belt amplitude time course
%                           for GE: RespData...
%
% OUT
%   cpulse              time events of R-wave peak in cardiac time series (seconds)
%                       for GE: usually empty
%   r                   respiratory time series
%   t                   vector of time points (in seconds)
%                       NOTE: This assumes the default sampling rate of 40
%                       Hz
%   c                   cardiac time series (ECG or pulse oximetry)
%
% EXAMPLE
%   [ons_secs.cpulse, ons_secs.rpulse, ons_secs.t, ons_secs.c] =
%       tapas_physio_read_physlogfiles_siemens(logfile, vendor, cardiac_modality);
%
%   See also tapas_physio_main_create_regressors
%
% Author: Lars Kasper
%         file structure information from PhLeM Toolbox, T. Verstynen (November 2007);
%                and Deshpande and J. Grinstead, Siemens Medical Solutions (March 2009)
%         additional log information Miriam Sebold, Charite Berlin (2014)
%
% Created: 2014-07-08
% Copyright (C) 2014 Institute for Biomedical Engineering, ETH/Uni Zurich.
%
% This file is part of the PhysIO toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id: tapas_physio_read_physlogfiles_siemens.m 466 2014-04-27 13:10:48Z kasperla $

%% read out values

if nargin < 2
    verbose.level = 0;
end
DEBUG = verbose.level >=2;

% process optional input parameters and overwrite defaults
defaults.ecgChannel         = 'mean'; % 'mean'; 'v1'; 'v2'
defaults.endCropSeconds     = 1;

args                = tapas_physio_propval(varargin, defaults);
tapas_physio_strip_fields(args);

cpulse              = [];
dt                  = log_files.sampling_interval;


if ~isempty(log_files.cardiac)
    
    [lineData, logFooter] = tapas_physio_read_physlogfiles_siemens_raw(...
        log_files.cardiac);
    
      
    % Determine relative start of acquisition from dicom headers and
    % logfile footers
    hasScanTimingDicomImage = ~isempty(log_files.scan_timing);
    
    if hasScanTimingDicomImage
          
        % load dicom
        useDicom = true;
        if useDicom
            dicomHeader             = spm_dicom_headers(...
                fullfile(log_files.scan_timing));
            
            ScanStartTimeSeconds    = dicomHeader{1}.AcquisitionTime;
            
            % TODO: Include AcquisitionNumber? InstanceNumber?
            ScanStopTimeSeconds     = dicomHeader{1}.AcquisitionTime + ...
                dicomHeader{1}.RepetitionTime/1000;
        else
            ScanStartTimeSeconds = logFooter.ScanStartTimeSeconds;
            ScanStopTimeSeconds = logFooter.ScanStopTimeSeconds;
        end
        
        switch log_files.align_scan
            case 'first'
                relative_start_acquisition = ScanStartTimeSeconds - ...
                    logFooter.LogStartTimeSeconds;
            case 'last'
                relative_start_acquisition = ScanStopTimeSeconds - ...
                    logFooter.LogStopTimeSeconds;
        end
    else
        relative_start_acquisition = 0;
    end           
        
    % add arbitrary offset specified by user
    relative_start_acquisition = relative_start_acquisition + ...
        log_files.relative_start_acquisition;
    
    data_table = tapas_physio_siemens_line2table(lineData);
    dataCardiac = tapas_physio_siemens_table2cardiac(data_table, ecgChannel, dt, ...
        relative_start_acquisition, endCropSeconds);
     
    if DEBUG
       verbose.fig_handles(end+1) = ...
           tapas_physio_plot_raw_physdata_siemens(dataCardiac);  
    end
    
    
    % crop end of log file
    cpulse = dataCardiac.cpulse_on;
    c = dataCardiac.c;
    t = dataCardiac.t;
    stopSample = dataCardiac.stopSample;
    cpulse(cpulse > t(dataCardiac.stopSample)) = [];
    t(stopSample+1:end) = [];
    c(stopSample+1:end) = [];
    
else
    c = [];
end

if ~isempty(log_files.respiration)
    r = load(log_files.respiration, 'ascii');
    nSamples = size(r,1);
    t = relative_start_acquisition + ((0:(nSamples-1))*dt)';
else
    r = [];
end

end
