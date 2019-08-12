function tests = tapas_physio_examples_test()
% Tests all PhysIO examples both as matlab scripts and via batch editor and
% compares relevant parts of physio-output structure and
% multiple_regressors file to reference results, saved by
% physio_update_examples when deploying
%
%    tests = tapas_physio_examples_test()
%
% NOTE: In order to run these tests, the corresponding example data
%       for PhysIO has be downloaded via tapas_download_example_data
%       (in misc/ subfolder of tapas release)
%
% IN
%
% OUT
%
% EXAMPLE
%   tapas_physio_examples_test
%
%   See also tapas_download_example_data physio_update_examples

% Author:   Lars Kasper
% Created:  2019-08-12
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

function test_ge_ppu3t(testCase)
%% Compares previously saved physio-structure and multiple regressors file
% to current output of re-run of GE PPU3T example
% Note: both SPM or matlab-script based execution is possible
% (check parameter doUseSpm below!)

% methods for recursively comparing structures, see
% https://ch.mathworks.com/help/matlab/ref/matlab.unittest.constraints.structcomparator-class.html
import matlab.unittest.constraints.IsEqualTo
import matlab.unittest.constraints.AbsoluteTolerance
import matlab.unittest.constraints.RelativeTolerance
import matlab.unittest.constraints.StructComparator
import matlab.unittest.constraints.NumericComparator
import matlab.unittest.constraints.StringComparator

doUseSpm = true;

% run GE example and extract physio
pathPhysioPublic = fullfile(fileparts(mfilename('fullpath')), '..', '..');
% TODO: Make generic!
pathExamples =  fullfile(pathPhysioPublic, '..', 'examples');

if doUseSpm
    pathCurrentExample = fullfile(pathExamples, 'GE/PPU3T');
    pathNow = pwd;
    cd(pathCurrentExample); % for prepending absolute paths correctly
    fileExample = fullfile(pathCurrentExample, 'ge_ppu3t_spm_job.mat');
    load(fileExample, 'matlabbatch');
    
    % remove unnecessary verbosity for test
    matlabbatch{1}.spm.tools.physio.verbose.level = 0;
    
    spm_jobman('run', matlabbatch);
    cd(pathNow);
    
    % retrieve physio struct from saved file
    matlabbatch{1}.spm.tools.physio.model.output_physio = fullfile(pathCurrentExample, ...
        matlabbatch{1}.spm.tools.physio.save_dir{1}, ...
        matlabbatch{1}.spm.tools.physio.model.output_physio);
    load(matlabbatch{1}.spm.tools.physio.model.output_physio, 'physio');
    actPhysio = physio;
else % has verbosity...cannot switch it off
    fileExample = fullfile(pathExamples, 'GE/PPU3T/ge_ppu3t_matlab_script.m');
    run(fileExample); % will output a PhysIO-struct
    actPhysio = physio;
end


% load physio from reference data
fileReferenceData = fullfile(pathExamples, 'TestReferenceResults', 'examples', ...
    'GE/PPU3T', 'physio.mat');
load(fileReferenceData, 'physio');
expPhysio = physio;

% extract cpulse from actual and expected solution and compare
actSolution = actPhysio.ons_secs.cpulse;
expSolution = expPhysio.ons_secs.cpulse;

verifyEqual(testCase, actSolution, expSolution);

%% compare all numeric sub-fields of physio with some tolerance

% ons_secs has all the computed preprocessed physiological and scan timing
% sync data, from which .model derives the physiological regressors later
% on
% spulse_per_vol cannot be compared, because cell!
testCase.verifyThat(actPhysio.ons_secs, ...
    IsEqualTo(expPhysio.ons_secs,  ...
    'Using', StructComparator(NumericComparator, 'Recursively', true), ...
    'Within', RelativeTolerance(0.01), ...
    'IgnoringFields',  {'spulse_per_vol'}...
    ));

% recursive with string
% testCase.verifyThat(actPhysio, ...
%     IsEqualTo(expPhysio, 'Using', ...
%     StructComparator(StringComparator, 'Recursively', true), ...
%     'IgnoringCase', true, ...
%     'IgnoringWhitespace', true ...
%     ));

end