function summary = runPublicationResults(opts)
% runPublicationResults  Run ManipulatorArm analyses needed for publication.
%
% Usage:
%   runPublicationResults()
%   runPublicationResults(struct('continueOnError', false))
%   runPublicationResults(struct('selectedTasks', ["testC", "testD"]))
%
% The child files are scripts and many start with "clear".  They are
% therefore executed in the base workspace so their clear statements do not
% destroy this orchestrator's function workspace.

    if nargin < 1 || isempty(opts)
        opts = struct();
    end
    opts = applyDefaultOptions(opts);

    scriptDir = fileparts(mfilename('fullpath'));
    projectRoot = fullfile(scriptDir, '..', '..');
    addpath(genpath(projectRoot));

    resultRoot = fullfile(scriptDir, 'results', 'publication_run');
    if ~exist(resultRoot, 'dir')
        mkdir(resultRoot);
    end

    runStamp = datestr(now, 'yyyymmdd_HHMMSS');
    logFile = fullfile(resultRoot, ['publication_run_' runStamp '.log']);
    summaryFile = fullfile(resultRoot, ['publication_summary_' runStamp '.csv']);
    latestSummaryFile = fullfile(resultRoot, 'publication_summary_latest.csv');

    tasks = publicationTaskList(scriptDir);
    tasks = filterTasks(tasks, opts);

    fprintf('Publication result orchestration\n');
    fprintf('  Result root      : %s\n', resultRoot);
    fprintf('  Log file         : %s\n', logFile);
    fprintf('  Tasks to execute : %d\n', numel(tasks));
    fprintf('  Continue on error: %d\n\n', opts.continueOnError);

    diary(logFile);
    diary on;

    rows = repmat(emptySummaryRow(), numel(tasks), 1);
    totalTic = tic;

    for k = 1:numel(tasks)
        task = tasks(k);
        fprintf('\n============================================================\n');
        fprintf('[%d/%d] %s\n', k, numel(tasks), char(task.label));
        fprintf('  Script : %s\n', char(task.script));
        fprintf('  Group  : %s\n', char(task.group));
        fprintf('============================================================\n');

        row = emptySummaryRow();
        row.order = k;
        row.id = string(task.id);
        row.label = string(task.label);
        row.group = string(task.group);
        row.script = string(task.script);
        row.resultDir = string(task.resultDir);
        row.status = "not_run";
        row.errorMessage = "";
        row.durationSec = NaN;
        row.expectedOutputsPresent = false;
        row.completedAt = "";

        taskTic = tic;
        try
            if opts.skipExisting && expectedOutputsPresent(task)
                fprintf('Skipping because expected outputs already exist.\n');
                row.status = "skipped_existing";
            else
                runScriptInBase(task.script);
                row.status = "ok";
            end
        catch ME
            row.status = "failed";
            row.errorMessage = string(ME.message);
            fprintf(2, '\nTask failed: %s\n', ME.message);
            fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
            if ~opts.continueOnError
                row.durationSec = toc(taskTic);
                row.expectedOutputsPresent = expectedOutputsPresent(task);
                row.completedAt = string(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
                rows(k) = row;
                writeSummary(rows(1:k), summaryFile, latestSummaryFile);
                diary off;
                rethrow(ME);
            end
        end

        row.durationSec = toc(taskTic);
        row.expectedOutputsPresent = expectedOutputsPresent(task);
        row.completedAt = string(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        rows(k) = row;
        writeSummary(rows(1:k), summaryFile, latestSummaryFile);

        fprintf('Task status: %s, duration %.1f s, expected outputs present: %d\n', ...
            char(row.status), row.durationSec, row.expectedOutputsPresent);
    end

    summary = struct2table(rows);
    fprintf('\n============================================================\n');
    fprintf('Publication orchestration complete in %.1f s\n', toc(totalTic));
    fprintf('  Summary: %s\n', summaryFile);
    fprintf('  Latest : %s\n', latestSummaryFile);
    fprintf('============================================================\n');

    compareStructuralMetrics(scriptDir, resultRoot);

    diary off;
end

function opts = applyDefaultOptions(opts)
    if ~isfield(opts, 'continueOnError'), opts.continueOnError = true; end
    if ~isfield(opts, 'skipExisting'), opts.skipExisting = true; end
    if ~isfield(opts, 'selectedTasks'), opts.selectedTasks = strings(0, 1); end
    if ~isfield(opts, 'selectedGroups'), opts.selectedGroups = strings(0, 1); end
end

function tasks = publicationTaskList(scriptDir)
    tasks = [
        % task('reference_loads', 'Reference load validation', 'validation', ...
        %     'testReferenceModuleLoads.m', '.', strings(0, 1))
        % 
        % task('reference_single_simp', 'Reference module single-load SIMP', 'reference_module', ...
        %     'runReferenceModuleSingleLoadSIMP.m', '.', ...
        %     ["referenceModuleSingleLoadSIMP.mat", "referenceModuleSingleLoad_summary.csv"])
        % 
        % task('reference_multi_simp', 'Reference module multi-load SIMP', 'reference_module', ...
        %     'runReferenceModuleMultiLoadSIMP.m', '.', ...
        %     ["referenceModuleMultiLoadSIMP.mat", "referenceModuleMultiLoad_summary.csv"])
        % 
        % task('reference_sweep_simp', 'Reference module SIMP sweep', 'reference_module', ...
        %     'sweepReferenceModuleSIMP.m', 'referenceModuleSIMP_sweep', ["summary.csv"])
        % 
        % task('stage2d_non_shell', 'Stage 2D non-shell topology', 'stage', ...
        %     'stage2D_nonShellTopology.m', 'stage2D_nonShellTopology', ["summary.csv"])
        % 
        % task('stage3a_verification', 'Stage 3A full-arm verification', 'stage', ...
        %     'stage3A_fullArmVerification.m', 'stage3A_fullArmVerification', ["summary.csv"])
        % 
        % task('stage3a_extended', 'Stage 3A extended metrics', 'stage', ...
        %     'stage3A_extendedMetrics.m', 'stage3A_fullArmVerification', ["extended_summary.csv"])
        % 
        % task('testA', 'Test A full-arm unlinked SIMP', 'full_arm_simp', ...
        %     'testA_fullArmUnlinkedSIMP.m', 'testA_fullArmUnlinkedSIMP', ["result.mat", "summary.csv"])

        task('testB', 'Test B full-arm linked SIMP', 'full_arm_simp', ...
            'testB_fullArmLinkedSIMP.m', 'testB_fullArmLinkedSIMP', ["result.mat", "summary.csv"])

        task('testE', 'Test E inverse full-arm unlinked SIMP', 'full_arm_simp_inverse', ...
            'testE_fullArmUnlinkedSIMP.m', 'testE_fullArmUnlinkedSIMP', ["result.mat", "summary.csv"])

        task('testF', 'Test F inverse full-arm linked SIMP', 'full_arm_simp_inverse', ...
            'testF_fullArmLinkedSIMP.m', 'testF_fullArmLinkedSIMP', ["result.mat", "summary.csv"])

        task('testG', 'Test G stress-constrained full-arm unlinked SIMP', 'full_arm_simp_stress', ...
            'testG_fullArmStressUnlinkedSIMP.m', 'testG_fullArmStressUnlinkedSIMP', ["result.mat", "summary.csv"])

        task('testH', 'Test H stress-constrained full-arm linked SIMP', 'full_arm_simp_stress', ...
            'testH_fullArmStressLinkedSIMP.m', 'testH_fullArmStressLinkedSIMP', ["result.mat", "summary.csv"])

        task('testC', 'Test C linked stress-intensity ESO', 'full_arm_stress_intensity', ...
            'testC_fullArmLinkedStressIntensity.m', 'testC_fullArmLinkedStressIntensity', ["result.mat", "history.csv"])

        task('testD', 'Test D unlinked stress-intensity ESO', 'full_arm_stress_intensity', ...
            'testD_fullArmUnlinkedStressIntensity.m', 'testD_fullArmUnlinkedStressIntensity', ["result.mat", "history.csv"])
    ];

    for k = 1:numel(tasks)
        tasks(k).script = fullfile(scriptDir, tasks(k).script);
        if tasks(k).resultDir == "."
            tasks(k).resultDir = fullfile(scriptDir, 'results');
        else
            tasks(k).resultDir = fullfile(scriptDir, 'results', tasks(k).resultDir);
        end
    end
end

function t = task(id, label, group, script, resultDir, expectedFiles)
    t.id = string(id);
    t.label = string(label);
    t.group = string(group);
    t.script = string(script);
    t.resultDir = string(resultDir);
    t.expectedFiles = string(expectedFiles(:));
end

function tasks = filterTasks(tasks, opts)
    selectedTasks = string(opts.selectedTasks(:));
    selectedGroups = string(opts.selectedGroups(:));

    if ~isempty(selectedTasks)
        keep = ismember([tasks.id], selectedTasks);
        tasks = tasks(keep);
    end

    if ~isempty(selectedGroups)
        keep = ismember([tasks.group], selectedGroups);
        tasks = tasks(keep);
    end

    if isempty(tasks)
        error('runPublicationResults:NoTasksSelected', ...
            'No publication tasks selected. Check selectedTasks/selectedGroups.');
    end
end

function runScriptInBase(scriptPath)
    scriptPath = char(scriptPath);
    escapedPath = strrep(scriptPath, '''', '''''');
    evalin('base', sprintf('run(''%s'');', escapedPath));
end

function ok = expectedOutputsPresent(task)
    expectedFiles = task.expectedFiles;
    expectedFiles(expectedFiles == "") = [];
    if isempty(expectedFiles)
        ok = false;
        return;
    end

    ok = exist(char(task.resultDir), 'dir') == 7;
    for k = 1:numel(expectedFiles)
        ok = ok && exist(char(fullfile(task.resultDir, expectedFiles(k))), 'file') == 2;
    end
end

function row = emptySummaryRow()
    row.order = NaN;
    row.id = "";
    row.label = "";
    row.group = "";
    row.script = "";
    row.resultDir = "";
    row.status = "";
    row.errorMessage = "";
    row.durationSec = NaN;
    row.expectedOutputsPresent = false;
    row.completedAt = "";
end

function writeSummary(rows, summaryFile, latestSummaryFile)
    summary = struct2table(rows);
    writetable(summary, summaryFile);
    writetable(summary, latestSummaryFile);
end

function compareStructuralMetrics(scriptDir, outputDir)
% Load structural_metrics.csv from each test result directory and print a
% cross-test comparison table.  Saves all_structural_metrics.csv.
    testIds = ["testA", "testB", "testC", "testD", "testG", "testH"];
    testDirs = [
        "testA_fullArmUnlinkedSIMP"
        "testB_fullArmLinkedSIMP"
        "testC_fullArmLinkedStressIntensity"
        "testD_fullArmUnlinkedStressIntensity"
        "testG_fullArmStressUnlinkedSIMP"
        "testH_fullArmStressLinkedSIMP"
    ];

    allRows = {};
    fprintf('\n============================================================\n');
    fprintf('Structural performance metrics comparison\n');
    fprintf('============================================================\n');
    fprintf('%-8s  %-14s  %10s  %10s  %10s  %10s  %8s\n', ...
        'Test', 'Config', 'sHM_ref', 'sHM_final', 'sHM_ratio', 'u_ratio', 'volFrac');
    fprintf('%s\n', repmat('-', 1, 75));

    for i = 1:numel(testIds)
        csvPath = fullfile(scriptDir, 'results', testDirs(i), 'structural_metrics.csv');
        if ~exist(char(csvPath), 'file')
            fprintf('%-8s  (no structural_metrics.csv found)\n', char(testIds(i)));
            continue;
        end
        T = readtable(char(csvPath));
        T.test_id = repmat(testIds(i), height(T), 1);
        allRows{end+1} = T;

        % Print ALL_MAX row for quick comparison
        maxRow = T(strcmp(T.config_name, 'ALL_MAX'), :);
        if isempty(maxRow), maxRow = T(end, :); end
        fprintf('%-8s  %-14s  %10.3e  %10.3e  %10.3f  %10.3f  %8.4f\n', ...
            char(testIds(i)), 'ALL_MAX', ...
            maxRow.sHM_ref(1), maxRow.sHM_final(1), ...
            maxRow.sHM_ratio(1), maxRow.u_ratio(1), maxRow.vol_frac(1));
    end

    if isempty(allRows)
        fprintf('No structural_metrics.csv files found. Run tests first.\n');
        return;
    end

    combined = vertcat(allRows{:});
    % Reorder columns: test_id first
    cols = combined.Properties.VariableNames;
    newOrder = ['test_id', cols(~strcmp(cols, 'test_id'))];
    combined = combined(:, newOrder);

    outFile = fullfile(outputDir, 'all_structural_metrics.csv');
    writetable(combined, outFile);
    fprintf('\nFull comparison saved to:\n  %s\n', outFile);
end
