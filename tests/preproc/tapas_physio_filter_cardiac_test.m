function tests = tapas_physio_filter_cardiac_test.m()
% Tests whether bandpass filter on PPU example data works as expected
%
%   tests = tapas_physio_filter_cardiac_test.m()
%
% IN
%
% OUT
%
% EXAMPLE
%   tapas_physio_filter_cardiac_test.m
%
%   See also
 
% Author:   Lars Kasper
% Created:  2019-07-02
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

function test_ge_ppu3t_filter(testCase)
%% Compares previously saved cpulse (detected cardiac pulses) from 
% physio-structure to same output when re-running current version of
% GE PPU3T example
% both SPM or matlab-script based execution is possible (check parameter
% doUseSpm below!)

% run GE example and extract physio
pathPhysioPublic = fullfile(fileparts(mfilename('fullpath')), '..', '..');
% TODO: Make generic!
pathExamples =  fullfile(pathPhysioPublic, '..', 'examples');

    pathCurrentExample = fullfile(pathExamples, 'Philips/PPU7T');
    pathNow = pwd;
    cd(pathCurrentExample); % for prepending absolute paths correctly
    fileExample = fullfile(pathCurrentExample, 'philips_ppu7t_spm_job.m');
    run(fileExample); % retrieve matlabbatch
    
    % remove unnecessary verbosity and processing of resp data
    matlabbatch{1}.spm.tools.physio.verbose.level = 0;
    matlabbatch{1}.spm.tools.physio.log_files.respiration = {''};

    spm_jobman('run', matlabbatch);
    cd(pathNow);
    
    % retrieve physio struct from saved file
    matlabbatch{1}.spm.tools.physio.model.output_physio = fullfile(pathCurrentExample, ...
        matlabbatch{1}.spm.tools.physio.save_dir{1}, ...
        matlabbatch{1}.spm.tools.physio.model.output_physio);
    load(matlabbatch{1}.spm.tools.physio.model.output_physio, 'physio');
    actPhysio = physio;
end
