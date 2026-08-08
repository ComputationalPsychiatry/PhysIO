classdef tapas_physio_readin_bids_missing_trigger_test < matlab.unittest.TestCase
    % Regression tests for BIDS recordings without an optional trigger column.

    methods (Test)
        function test_unified_file_without_trigger(testCase)
            pathTestData = tempname;
            mkdir(pathTestData);
            testCase.addTeardown(@() rmdir(pathTestData, 's'));

            fileTsv = fullfile(pathTestData, 'minimal_physio.tsv');
            fid = fopen(fileTsv, 'w');
            fprintf(fid, '1\t2\n3\t4\n5\t6\n');
            fclose(fid);

            fileJson = fullfile(pathTestData, 'minimal_physio.json');
            fid = fopen(fileJson, 'w');
            fprintf(fid, ['{"SamplingFrequency":1000,"StartTime":0,' ...
                '"Columns":["cardiac","respiratory"]}']);
            fclose(fid);

            physio = tapas_physio_new();
            physio.log_files.cardiac = fileTsv;
            physio.log_files.respiration = fileTsv;
            physio.log_files.scan_timing = fileJson;
            physio.log_files.sampling_interval = [];
            physio.log_files.relative_start_acquisition = [];
            physio.verbose.level = 0;

            [c, r, t, cpulse, acq_codes] = ...
                tapas_physio_read_physlogfiles_bids(physio.log_files, ...
                'PPU', physio.verbose);

            testCase.verifyEqual(c, [1; 3; 5]);
            testCase.verifyEqual(r, [2; 4; 6]);
            testCase.verifyEqual(t, [0; 0.001; 0.002], 'AbsTol', eps);
            testCase.verifyEmpty(cpulse);
            testCase.verifyEmpty(acq_codes);
        end
    end
end
