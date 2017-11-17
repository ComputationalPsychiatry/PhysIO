function [cTics, c, cSignals, extTriggerSignals, stringChannels] = ...
                tapas_physio_split_data_per_channel_siemens_tics(C, ecgChannel)
% Splits column data of Siemens Tics format according to their channel
% label in the 2nd column
%
%      [cTics, c, cSignals, extTriggerSignals, stringChannels] = ...
%                tapas_physio_split_data_per_channel_siemens_tics(C, ecgChannel);
%        
%
%
%   Example structure of a tics logfile (CMRR):
% 
%      19824786     ECG2   2084 
%      19824787     ECG2   2190 
%      19824788     ECG2   2198  PULS_TRIGGER
%      19824789     ECG2   2095 
%      ...
%      19824762     ECG3   1948 
%      19824763     ECG3   1940 
%      19824764     ECG3   1953 
%      
% IN
%
% OUT
%
% EXAMPLE
%   tapas_physio_split_data_per_channel_siemens_tics
%
%   See also
%
% Author: Lars Kasper
% Created: 2017-11-17
% Copyright (C) 2017 TNU, Institute for Biomedical Engineering,
%                    University of Zurich and ETH Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under the terms of the GNU General Public
% License (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id: teditRETRO.m 775 2015-07-17 10:52:58Z kasperla $


  uniqueCTics = unique(cTics);
            
            %% retrieve data channel-wise, and interpolate to common grid
            for iChannel = 1:nChannels
                idxChannel = find(strcmp(C{2}, stringChannels{iChannel}));
                cChannel = c(idxChannel);
                cTicsChannel = cTics(idxChannel);
                
                %% interpolate to same time grid (tics) for channel combination  already...
                
                % first, remove duplicates in tics time axis by averaging
                % their corresponding values, keep only one occurence
                
                % detects bins with more than one entry, index of first
                % occurence is returned
                idxDuplicateTics = find(hist(cTicsChannel,unique(cTicsChannel))>1);
                
                idxDuplicatesAll = [];
                for idx = idxDuplicateTics
                    % detect all occurences that match first labeled as
                    % duplicate
                    idxDuplicates = reshape(find(cTicsChannel==cTicsChannel(idx)), 1, []);
                    cChannel(idxDuplicates(1)) = mean(cChannel(idxDuplicateTics));
                    
                    % label for later deletion without changing index order
                    % just now
                   idxDuplicatesAll = [idxDuplicatesAll, idxDuplicates(2:end)];
                end
                 
                % remove duplicates
                cChannel(idxDuplicatesAll) = [];
                cTicsChannel(idxDuplicatesAll) = [];
                idxChannel(idxDuplicatesAll) = [];
                
                % now interpolate without duplicates!
                cChannelInterpolated = interp1(cTicsChannel, cChannel, ...
                    uniqueCTics, 'linear', 'extrap');
                
                % save interpolated data
                indPerChannel{iChannel} = idxChannel;
                cPerChannel{iChannel} = cChannel;
                cTicsPerChannel{iChannel} = cTicsChannel;
                cPerChannelInterpolated{iChannel} = cChannelInterpolated;

            end
            
            switch ecgChannel
                case 'v1'
                    
                case 'v2'
                case 'mean'
            end