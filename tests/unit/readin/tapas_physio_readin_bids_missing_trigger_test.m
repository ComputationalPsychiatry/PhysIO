classdef tapas_physio_readin_bids_missing_trigger_test < matlab.unittest.TestCase
    % Regression tests for BIDS recordings without an optional trigger column.

    methods (Test)
        function test_unified_file_without_trigger(testCase)
            [logFiles, verbose] = testCase.createBidsFixture( ...
                '1\t2\n3\t4\n5\t6\n', ...
                {'cardiac', 'respiratory'});

            [c, r, t, cpulse, acq_codes] = ...
                tapas_physio_read_physlogfiles_bids(logFiles, ...
                'PPU', verbose);

            testCase.verifyEqual(c, [1; 3; 5]);
            testCase.verifyEqual(r, [2; 4; 6]);
            testCase.verifyEqual(t, [0; 0.001; 0.002], 'AbsTol', eps);
            testCase.verifyEmpty(cpulse);
            testCase.verifyEmpty(acq_codes);
        end

        function test_unified_file_with_trigger(testCase)
            [logFiles, verbose] = testCase.createBidsFixture( ...
                '1\t2\t1\n3\t4\t0\n5\t6\t1\n', ...
                {'cardiac', 'respiratory', 'trigger'});

            [c, r, t, cpulse, acq_codes] = ...
                tapas_physio_read_physlogfiles_bids(logFiles, ...
                'PPU', verbose);

            testCase.verifyEqual(c, [1; 3; 5]);
            testCase.verifyEqual(r, [2; 4; 6]);
            testCase.verifyEqual(t, [0; 0.001; 0.002], 'AbsTol', eps);
            testCase.verifyEmpty(cpulse);
            testCase.verifyEqual(acq_codes, [8; 0; 8]);
        end
    end

    methods (Access = private)
        function [logFiles, verbose] = createBidsFixture( ...
                testCase, tsvContents, columns)
            pathTestData = tempname;
            mkdir(pathTestData);
            testCase.addTeardown(@() rmdir(pathTestData, 's'));

            fileTsv = fullfile(pathTestData, 'minimal_physio.tsv');
            fid = fopen(fileTsv, 'w');
            fprintf(fid, tsvContents);
            fclose(fid);

            fileJson = fullfile(pathTestData, 'minimal_physio.json');
            metadata = struct('SamplingFrequency', 1000, ...
                'StartTime', 0, 'Columns', {columns});
            fid = fopen(fileJson, 'w');
            fprintf(fid, '%s', jsonencode(metadata));
            fclose(fid);

            physio = tapas_physio_new();
            physio.log_files.cardiac = fileTsv;
            physio.log_files.respiration = fileTsv;
            physio.log_files.scan_timing = fileJson;
            physio.log_files.sampling_interval = [];
            physio.log_files.relative_start_acquisition = [];
            physio.verbose.level = 0;

            logFiles = physio.log_files;
            verbose = physio.verbose;
        end
    end
end
