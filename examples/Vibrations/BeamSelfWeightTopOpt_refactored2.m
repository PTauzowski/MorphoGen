%% Beam Self-Weight Topology Optimization
% Main script for topology optimization of beams under self-weight and 
% harmonic loading conditions with natural frequency analysis

clearvars;
close all;
clc;

%% Configuration
CONFIG = struct();
CONFIG.dataFile = 'BeamSelfWeightTopOpt.json';
CONFIG.resultsFile = 'BeamSelfWeightTopOptConst.mat';
CONFIG.multiFormsResultsFile = 'BeamMultiFormsTopOptConst.mat';
CONFIG.resultsFolder = 'beam_results_const_ref';
CONFIG.forceRecompute = false;
CONFIG.nEigenForms = 30;
CONFIG.nFormsToPlot = 10;
CONFIG.loadModesNumber = 3;
CONFIG.alphaDivisions = 4;

%% Initialize
data = jsondecode(fileread(CONFIG.dataFile));
createDirectoryIfNeeded(CONFIG.resultsFolder);

%% Part 1: Self-Weight Topology Optimization
if CONFIG.forceRecompute || ~exist(CONFIG.resultsFile, "file")
    fprintf('Computing self-weight topology optimization...\n');
    [topOpt, mesh, analysisLinear, constElems] = computeSelfWeightTopology(data, CONFIG);
    % Save only computation results, not CONFIG
    save(CONFIG.resultsFile, 'topOpt', 'mesh', 'analysisLinear', 'constElems', 'data');
    fprintf('Self-weight optimization complete.\n');
else
    fprintf('Loading existing self-weight results...\n');
    % Load only computation results, CONFIG remains from current script
    load(CONFIG.resultsFile, 'topOpt', 'mesh', 'analysisLinear', 'constElems', 'data');
end

%% Part 2: Natural Vibration Analysis for Design Domain
fprintf('Analyzing design domain natural vibrations...\n');
[dd_freq, dd_forms] = analyzeDesignDomainVibrations(analysisLinear, mesh, CONFIG);
plotDesignDomainResults(dd_freq, dd_forms, mesh, analysisLinear, CONFIG);

%% Part 3: Self-Weight Topology Natural Vibrations
fprintf('Analyzing self-weight topology vibrations...\n');
analyzeSelfWeightTopology(topOpt, analysisLinear, mesh, data, dd_forms, CONFIG);

%% Part 4: Multi-Form Topology Optimization
if CONFIG.forceRecompute || ~exist(CONFIG.multiFormsResultsFile, "file")
    fprintf('Computing multi-form topology optimizations...\n');
    results = computeMultiFormTopologies(data, mesh, analysisLinear, constElems, CONFIG);
    % Save only computation results, not CONFIG
    save(CONFIG.multiFormsResultsFile, '-struct', 'results', '-v7.3');
    fprintf('Multi-form optimization complete.\n');
else
    fprintf('Loading existing multi-form results...\n');
    % Load computation results, CONFIG remains from current script
    results = load(CONFIG.multiFormsResultsFile);
end

%% Part 5: Generate All Output Results
fprintf('Generating output visualizations...\n');
generateAllOutputs(results, dd_forms, data, CONFIG);
fprintf('Analysis complete!\n');

%% ========================================================================
%% FUNCTION DEFINITIONS
%% ========================================================================

function [topOpt, mesh, analysisLinear, constElems] = computeSelfWeightTopology(data, CONFIG)
    % Setup mesh and finite elements
    [mesh, fe, material] = setupMeshAndElements(data);
    
    % Setup linear analysis with boundary conditions
    analysisLinear = setupLinearAnalysis(fe, mesh, data);
    analysisLinear.selfLoadFactor = 1;
    
    % Define passive elements (flanges)
    constElems = definePassiveElements(mesh, data);
    
    % Create and solve topology optimization
    topOpt = createTopologyOptimization(data, analysisLinear, constElems);
    [~, ~] = topOpt.solve();
    close all;
end

function [mesh, fe, material] = setupMeshAndElements(data)
    % Create mesh
    sfL4 = ShapeFunctionQ4;
    mesh = Mesh();
    
    h = data.domain.size.height;
    l = data.domain.size.length;
    nelx = data.domain.mesh.nelx;
    nely = data.domain.mesh.nely;
    
    mesh.addRectMesh2D(0, 0, l, h, nelx, nely, sfL4.pattern);
    
    % Create finite element
    fe = PlaneStressElem(sfL4, mesh.elems);
    
    % Create and assign material
    material = PlaneStressMaterial('mat1');
    material.setElasticIzo(data.materials.planeStressIsotropic.E, ...
                          data.materials.planeStressIsotropic.nu);
    material.setMassIzoMatrix(data.materials.planeStressIsotropic.rho);
    fe.setMaterial(material);
end

function analysisLinear = setupLinearAnalysis(fe, mesh, data)
    analysisLinear = LinearElasticity(fe, mesh, false);
    
    l = data.domain.size.length;
    
    % Define boundary conditions
    leftEdge = Selector(@(x) abs(x(:,1)) < 0.0005);
    rightEdge = Selector(@(x) abs(x(:,1) - l) < 0.0005);
    
    analysisLinear.fixNodes(leftEdge, ["ux", "uy"]);
    analysisLinear.fixNodes(rightEdge, ["ux", "uy"]);
end

function constElems = definePassiveElements(mesh, data)
    res = data.domain.mesh.nely;
    l = data.domain.size.length;
    constThickness = round(res * data.domain.passive_regions.flange.top_relative_height);
    
    constRows = constThickness;
    ncel = round(constRows * res * l);
    nelems = size(mesh.elems, 1);
    
    constElems = [1:ncel, nelems:-1:(nelems - ncel)];
end

function topOpt = createTopologyOptimization(data, analysis, constElems)
    Rfilter = data.optimisation.filter_radius;
    cutThreshold = data.optimisation.cut_threshold;
    penal = data.optimisation.penalization;
    volFrac = min(data.optimisation.volume_fraction);
    
    topOpt = StressIntensityTopologyOptimizationVol(...
        Rfilter, analysis, cutThreshold, penal, volFrac, true);
    topOpt.setConstElems(constElems);
end

function [dd_freq, dd_forms] = analyzeDesignDomainVibrations(analysisLinear, mesh, CONFIG)
    vibrations = LinearNaturalVibration(analysisLinear.felems, mesh);
    vibrations.supports = analysisLinear.supports;
    vibrations.solve(CONFIG.nEigenForms, 1);
    
    dd_freq = vibrations.frequencies;
    dd_forms = vibrations.qforms;
end

function plotDesignDomainResults(dd_freq, dd_forms, mesh, analysisLinear, CONFIG)
    folderName = fullfile(CONFIG.resultsFolder, 'design_domain');
    createDirectoryIfNeeded(folderName);
    
    nelems = size(mesh.elems, 1);
    vibrations = LinearNaturalVibration(analysisLinear.felems, mesh);
    vibrations.supports = analysisLinear.supports;
    vibrations.frequencies = dd_freq;
    vibrations.qforms = dd_forms;
    
    filename = fullfile(folderName, 'design_domain');
    vibrations.plotNaturalForms(filename, 1:CONFIG.nFormsToPlot, 1:nelems, ...
                                '', 0.05, 'k', [0.8 0.8 0.8]);
end

function analyzeSelfWeightTopology(topOpt, analysisLinear, mesh, data, dd_forms, CONFIG)
    mainFolder = fullfile(CONFIG.resultsFolder, 'self_weight');
    createDirectoryIfNeeded(mainFolder);
    
    volumeFractions = data.optimisation.volume_fraction;
    topOpt.setFrame(topOpt.findFrame(volumeFractions));
    
    % Solve vibrations for optimized topology
    vibrations = LinearNaturalVibration(analysisLinear.felems, mesh);
    vibrations.supports = analysisLinear.supports;
    vibrations.solve(CONFIG.nEigenForms, topOpt.x);
    
    % Plot topology
    figure;
    topOpt.plotCurrentFrame();
    savefig(gcf, fullfile(mainFolder, 'self_weight_optimal_topology.fig'));
    
    % Plot natural forms
    vibrations.plotNaturalForms(...
        fullfile(mainFolder, 'self_weight_optimal_topology'), ...
        1:10, topOpt.x > 0.5, 'Self weight topology, ', 0.02, 'k', 'k');
    
    % Compute and save correlation matrix
    correlMatrix = vibrations.ComputeCorrelationMatrix(12, dd_forms(:, 1:3))';
    vibrations.printCorrelationTable(...
        fullfile(mainFolder, 'self_weight_topology_correlation_table'), ...
        'self weight topology correlation', correlMatrix);
end

function results = computeMultiFormTopologies(data, mesh, analysisLinear, constElems, CONFIG)
    results = struct();
    
    % Single mode optimizations
    fprintf('  Computing single mode topologies...\n');
    [results.formTopOpts, results.formVibrations] = ...
        computeSingleModeTopologies(data, mesh, analysisLinear, constElems, CONFIG);
    
    % Mixed mode optimizations
    alphas = (1:(CONFIG.alphaDivisions-1)) / CONFIG.alphaDivisions;
    results.alphas = alphas;
    
    fprintf('  Computing mixed mode topologies (1-2 envelope)...\n');
    [results.mixedTopOpts12env, results.mixedVibrations12env] = ...
        computeMixedModeTopologies(data, mesh, analysisLinear, constElems, ...
                                   [1, 2], alphas, 'envelope', CONFIG);
    
    fprintf('  Computing mixed mode topologies (1-2 average)...\n');
    [results.mixedTopOpts12av, results.mixedVibrations12av] = ...
        computeMixedModeTopologies(data, mesh, analysisLinear, constElems, ...
                                   [1, 2], alphas, 'average', CONFIG);
    
    fprintf('  Computing mixed mode topologies (1-3 envelope)...\n');
    [results.mixedTopOpts13env, results.mixedVibrations13env] = ...
        computeMixedModeTopologies(data, mesh, analysisLinear, constElems, ...
                                   [1, 3], alphas, 'envelope', CONFIG);
    
    fprintf('  Computing mixed mode topologies (1-3 average)...\n');
    [results.mixedTopOpts13av, results.mixedVibrations13av] = ...
        computeMixedModeTopologies(data, mesh, analysisLinear, constElems, ...
                                   [1, 3], alphas, 'average', CONFIG);
    
    close all;
end

function [topOpts, vibrations] = computeSingleModeTopologies(data, mesh, analysisLinear, constElems, CONFIG)
    topOpts = [];
    vibrations = [];
    
    Rfilter = data.optimisation.filter_radius;
    cutThreshold = data.optimisation.cut_threshold;
    penal = data.optimisation.penalization;
    volFrac = data.optimisation.volume_fraction;
    
    fe = analysisLinear.felems;
    
    for k = 1:CONFIG.loadModesNumber
        analysisHarmonic = ElasticHarmonicVibrations(fe, mesh, k, false, true);
        analysisHarmonic.supports = analysisLinear.supports;
        
        vib = LinearNaturalVibration(analysisLinear.felems, mesh);
        vib.supports = analysisLinear.supports;
        
        topOpt = StressIntensityTopologyOptimizationVol(...
            Rfilter, analysisHarmonic, cutThreshold, penal, volFrac, true);
        topOpt.setConstElems(constElems);
        
        figure; hold on;
        [~, xopt] = topOpt.solve();
        vib.solve(CONFIG.nEigenForms, xopt);
        
        vibrations = [vibrations, vib];
        topOpts = [topOpts, topOpt];
    end
end

function [topOpts, vibrations] = computeMixedModeTopologies(data, mesh, analysisLinear, constElems, modes, alphas, method, CONFIG)
    topOpts = [];
    vibrations = [];
    
    fe = analysisLinear.felems;
    
    % Create harmonic analyses for specified modes
    harmonicAnalyses = cell(length(modes), 1);
    for i = 1:length(modes)
        harmonicAnalyses{i} = ElasticHarmonicVibrations(fe, mesh, modes(i), false, true);
        harmonicAnalyses{i}.supports = analysisLinear.supports;
    end
    
    Rfilter = data.optimisation.filter_radius;
    cutThreshold = data.optimisation.cut_threshold;
    penal = data.optimisation.penalization;
    volFrac = data.optimisation.volume_fraction;
    
    for k = 1:length(alphas)
        vib = LinearNaturalVibration(analysisLinear.felems, mesh);
        vib.supports = analysisLinear.supports;
        
        weights = [1 - alphas(k), alphas(k)];
        
        if strcmp(method, 'envelope')
            topOpt = StressIntensityMultiMaxTopologyOptimization(...
                Rfilter, [harmonicAnalyses{:}], weights, cutThreshold, penal, volFrac, true);
        else % average
            topOpt = StressIntensityMultiAvTopologyOptimization(...
                Rfilter, [harmonicAnalyses{:}], weights, cutThreshold, penal, volFrac, true);
        end
        
        topOpt.setConstElems(constElems);
        
        figure; hold on;
        [~, xopt] = topOpt.solve();
        vib.solve(CONFIG.nEigenForms, xopt);
        
        topOpts = [topOpts, topOpt];
        vibrations = [vibrations, vib];
    end
end

function generateAllOutputs(results, dd_forms, data, CONFIG)
    volFrac = data.optimisation.volume_fraction;
    
    % Single mode results
    for formIdx = 1:CONFIG.loadModesNumber
        generateSingleModeOutput(results.formTopOpts(formIdx), ...
                                results.formVibrations(formIdx), ...
                                formIdx, volFrac, dd_forms, CONFIG);
    end
    
    % Mixed mode results
    generateMixedModeOutputs(results.mixedTopOpts12env, results.mixedVibrations12env, ...
                            [1, 2], 'envelope', results.alphas, volFrac, dd_forms, CONFIG);
    generateMixedModeOutputs(results.mixedTopOpts12av, results.mixedVibrations12av, ...
                            [1, 2], 'average', results.alphas, volFrac, dd_forms, CONFIG);
    generateMixedModeOutputs(results.mixedTopOpts13env, results.mixedVibrations13env, ...
                            [1, 3], 'envelope', results.alphas, volFrac, dd_forms, CONFIG);
    generateMixedModeOutputs(results.mixedTopOpts13av, results.mixedVibrations13av, ...
                            [1, 3], 'average', results.alphas, volFrac, dd_forms, CONFIG);
    
    % Generate combined visualizations
    fprintf('Generating combined visualizations...\n');
    generateCombinedVisualization(results.formTopOpts, results.mixedTopOpts12env, ...
                                 [1, 2], 'envelope', results.alphas, volFrac, CONFIG);
    generateCombinedVisualization(results.formTopOpts, results.mixedTopOpts12av, ...
                                 [1, 2], 'average', results.alphas, volFrac, CONFIG);
    generateCombinedVisualization(results.formTopOpts, results.mixedTopOpts13env, ...
                                 [1, 3], 'envelope', results.alphas, volFrac, CONFIG);
    generateCombinedVisualization(results.formTopOpts, results.mixedTopOpts13av, ...
                                 [1, 3], 'average', results.alphas, volFrac, CONFIG);
    
    close all;
end

function generateSingleModeOutput(topOpt, vibration, formIdx, volFrac, dd_forms, CONFIG)
    figure; hold on;
    iters = topOpt.findFrame(volFrac);
    topOpt.setFrame(iters);
    topOpt.plotCurrentFrame();
    title(sprintf('FSD, topology for Φ_%d, vol = %.2f, iteration = %d', ...
                  formIdx, volFrac(1), iters));
    
    folderName = fullfile(CONFIG.resultsFolder, sprintf('topology_for_%d_mode', formIdx));
    createDirectoryIfNeeded(folderName);
    
    savefig(gcf, fullfile(folderName, sprintf('topology_for_%d_mode.fig', formIdx)));
    
    vibration.plotNaturalForms(...
        fullfile(folderName, sprintf('topology_for_%d_mode', formIdx)), ...
        1:12, topOpt.x > 0.5, sprintf(', vol = %.2f, ', volFrac(1)), ...
        0.02, 'k', 'k');
    
    correlMatrix = vibration.ComputeCorrelationMatrix(12, dd_forms(:, 1:3))';
    vibration.printCorrelationTable(...
        fullfile(folderName, sprintf('correlations_for_%d_mode', formIdx)), ...
        'topology correlation', correlMatrix);
end

function generateMixedModeOutputs(topOpts, vibrations, modes, method, alphas, volFrac, dd_forms, CONFIG)
    modeStr = sprintf('%d_%d', modes(1), modes(2));
    folderName = fullfile(CONFIG.resultsFolder, sprintf('topology_for_modes_%s', modeStr));
    createDirectoryIfNeeded(folderName);
    
    for k = 1:length(alphas)
        alpha = alphas(k);
        
        figure; hold on;
        topOpts(k).setFrame(topOpts(k).findFrame(volFrac));
        topOpts(k).plotCurrentFrame();
        title(sprintf('FSD P=[%.2f*Φ_%d %.2f*Φ_%d]', ...
                     1-alpha, modes(1), alpha, modes(2)));
        
        filePrefix = sprintf('topology_for_modes_%s_%s', modeStr, method);
        subFolder = fullfile(CONFIG.resultsFolder, sprintf('%s_%d', filePrefix, k));
        createDirectoryIfNeeded(subFolder);
        
        savefig(gcf, fullfile(folderName, sprintf('%s_%d.fig', filePrefix, k)));
        
        vibrations(k).plotNaturalForms(...
            fullfile(subFolder, sprintf('%s_%d', filePrefix, k)), ...
            1:10, topOpts(k).x > 0.5, ...
            sprintf(' %s, vol=%.2f [%.2f*Φ_%d %.2f*Φ_%d] ', ...
                   method, volFrac(1), 1-alpha, modes(1), alpha, modes(2)), ...
            0.02, 'k', 'k');
        
        correlMatrix = vibrations(k).ComputeCorrelationMatrix(12, dd_forms(:, 1:3))';
        vibrations(k).printCorrelationTable(...
            fullfile(subFolder, sprintf('%s_%d', filePrefix, k)), ...
            'topology correlation', correlMatrix);
    end
end

function createDirectoryIfNeeded(dirName)
    if ~exist(dirName, 'dir')
        mkdir(dirName);
    end
end

function generateCombinedVisualization(formTopOpts, mixedTopOpts, modes, method, alphas, volFrac, CONFIG)
    % Create a figure with 5 subplots (5 rows, 1 column) showing: Phi1, mixed1, mixed2, mixed3, Phi2 (or Phi3)
    % For modes [1,2]: shows Phi1 -> combinations -> Phi2
    % For modes [1,3]: shows Phi1 -> combinations -> Phi3
    
    % Calculate number of plots (2 pure modes + number of mixed)
    nMixed = length(alphas);
    nTotal = 2 + nMixed;
    
    % Create main figure
    fig = figure('Position', [100, 100, 600, 1400]);
    
    % Plot first pure mode (Phi_modes(1))
    subplot(nTotal, 1, 1);
    formTopOpts(modes(1)).setFrame(formTopOpts(modes(1)).findFrame(volFrac));
    % Get current axes before plotting
    ax1 = gca;
    formTopOpts(modes(1)).plotCurrentFrame();
    % Restore this subplot as current
    subplot(nTotal, 1, 1);
    title(sprintf('\\Phi_{%d}', modes(1)), 'FontSize', 12, 'FontWeight', 'bold');
    axis equal tight;
    colormap(gca, 'gray');
    
    % Plot mixed modes
    for k = 1:nMixed
        subplot(nTotal, 1, k + 1);
        mixedTopOpts(k).setFrame(mixedTopOpts(k).findFrame(volFrac));
        mixedTopOpts(k).plotCurrentFrame();
        % Restore this subplot as current
        subplot(nTotal, 1, k + 1);
        alpha = alphas(k);
        title(sprintf('%.2f\\Phi_{%d} + %.2f\\Phi_{%d}', ...
                     1-alpha, modes(1), alpha, modes(2)), ...
              'FontSize', 12, 'FontWeight', 'bold');
        axis equal tight;
        colormap(gca, 'gray');
    end
    
    % Plot second pure mode (Phi_modes(2))
    subplot(nTotal, 1, nTotal);
    formTopOpts(modes(2)).setFrame(formTopOpts(modes(2)).findFrame(volFrac));
    formTopOpts(modes(2)).plotCurrentFrame();
    % Restore this subplot as current
    subplot(nTotal, 1, nTotal);
    title(sprintf('\\Phi_{%d}', modes(2)), 'FontSize', 12, 'FontWeight', 'bold');
    axis equal tight;
    colormap(gca, 'gray');
    
    % Add overall title
    sgtitle(sprintf('Topology Optimization: \\Phi_{%d} to \\Phi_{%d} (%s, vol=%.2f)', ...
                   modes(1), modes(2), method, volFrac(1)), ...
           'FontSize', 14, 'FontWeight', 'bold');
    
    % Adjust spacing between subplots
    set(fig, 'Units', 'normalized');
    
    % Save the combined figure
    modeStr = sprintf('%d_%d', modes(1), modes(2));
    folderName = fullfile(CONFIG.resultsFolder, sprintf('topology_for_modes_%s', modeStr));
    createDirectoryIfNeeded(folderName);
    
    filename = fullfile(folderName, sprintf('combined_%s_modes_%s.fig', method, modeStr));
    savefig(fig, filename);
    
    % Also save as PNG for easier viewing
    filename_png = fullfile(folderName, sprintf('combined_%s_modes_%s.png', method, modeStr));
    print(fig, filename_png, '-dpng', '-r300');
    
    fprintf('  Saved combined visualization: %s\n', filename_png);
end