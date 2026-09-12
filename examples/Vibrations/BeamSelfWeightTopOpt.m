clear;
clear classes;
clear all;
rehash
close all;


data = jsondecode(fileread('BeamSelfWeightTopOpt.json'));

baseDir = fileparts(mfilename('fullpath'));
%results_folder_name = fullfile(baseDir, "beam_results_const");
results_folder_name = "beam_results_const";

if ~isfolder(results_folder_name)
    mkdir(results_folder_name);
end

results_filename = 'BeamSelfWeightTopOptConst.mat';

% If true: always compute topology (ignore existing results)
% If false: compute only if the results file defined in results_filename does not exist
enforce_topology_computation = false;

if enforce_topology_computation || ~exist(results_filename, "file")

    % Cantilever topology optimization elastic task.
    
    % Resolution of shortest (vertical) edge
    res = data.domain.mesh.nely;
    
    % height of the cantilever
    h = data.domain.size.height;
    
    
    % length of the cantilever
    l = data.domain.size.length;
    
    % Filtering radius
    Rfilter = data.optimisation.filter_radius; % 1.5*h/res;
    
    %Removal intensity threshold
    cutTreshold = data.optimisation.cut_threshold; 
    
    %penalty factor
    penal = data.optimisation.penalization;
    
    % Type of shape function to be used (here: four node Langrange)
    sfL4 = ShapeFunctionQ4;
    
    % Creating FE mesh object
    mesh = Mesh();
    
    % Generating rectangular mesh ( aspect*h x h )
    mesh.addRectMesh2D(0, 0, l, h, data.domain.mesh.nelx, res, sfL4.pattern);
    nelems=size(mesh.elems,1);
    
    
    % Create plane stress finite element object
    fe = PlaneStressElem( sfL4, mesh.elems );
    
    E = data.materials.planeStressIsotropic.E;
    nu = data.materials.planeStressIsotropic.nu;
    rho = data.materials.planeStressIsotropic.rho;
    
    % Create isotropic material object
    material = PlaneStressMaterial('mat1');
    material.setElasticIzo(E, nu);
    material.setMassIzoMatrix(rho);
    % Assigning material to finite element
    fe.setMaterial( material );
    
    % Creating linear elastic finite element analysis object with weighted matrix feature, weighted by element density.
    analysisLinear = LinearElasticity( fe, mesh, false );
    
    % Creating node selector object to select fixed edge (left)
    %loadEdgeSelector = Selector( @(x)( abs(x(:,2) - h ) < 0.0005 ) );
    LeftEdgeSelector = Selector( @(x)( abs(x(:,1)) < 0.0005 ) );
    RightEdgeSelector = Selector( @(x)( abs(x(:,1)-l) < 0.0005 ) );
    
    % Fixing edges structure according to above defined node selector object
    analysisLinear.fixNodes( LeftEdgeSelector, ["ux" "uy"] );
    analysisLinear.fixNodes( RightEdgeSelector, ["ux" "uy"] );
    
    % Passive elements (height of flange)
    const_thickness = round(res*data.domain.passive_regions.flange.top_relative_height); 

    volumeFractions = data.optimisation.volume_fraction;

    const_rows = const_thickness;
    ncel = round(const_rows*res*l);
    const_elems = [1:ncel size(mesh.elems,1):-1:size(mesh.elems,1)-ncel ];
    
    % set load factor to enable self weight computation
    analysisLinear.selfLoadFactor=1;
            
    topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysisLinear, cutTreshold, penal, min( volumeFractions ), true );
    topOpt.setConstElems(const_elems);


    % Enable optional Woodbury reanalysis:
    % args: (useFlag, rankMax, eigTol, fallbackMode, verbose, monitorPrint)
    % topOpt.configureWoodburyReanalysis(true, 300, 1e-10, 'direct', true, true);
    % 
    % % Optional tuning:
    % topOpt.woodburyRebaseOnFallback = true;  % reset base matrix after fallback
    % topOpt.woodburyFallbackMode = 'pcg';    % optional PCG fallback
    
    
    % Uncomment if topology result have to be computed
    tic;
    [objF, xopt]  = topOpt.solve();
    toc
    %topOpt.printWoodburyMonitorSummary(20);

    close all;
    save(results_filename);

else

    load(results_filename);
    %figure;
    %topOpt.plotCurrentFrame();

end

results_folder_name = "beam_results_const";

%natural vibration problem setup
nEigenForms=30;
nforms=10;

% natural vibration analysis object definition
vibrations = LinearNaturalVibration( analysisLinear.felems, mesh );
vibrations.supports = analysisLinear.supports;

% solve natural vibration problem for design domain beam
vibrations.solve( nEigenForms, 1 ); 
dd_freq = vibrations.frequencies;
dd_forms = vibrations.qforms;

folderName = fullfile(results_folder_name,"design_domain");

if ~isfolder(folderName)
        mkdir(folderName);
end

filename_predictor = fullfile( folderName, "design_domain");
vibrations.plotNaturalForms( filename_predictor, 1:nforms,1:nelems,"",0.05,'k',[0.8 0.8 0.8]);


mainFolderName = fullfile(results_folder_name,"self_weight");
if ~isfolder( mainFolderName)
    mkdir( mainFolderName);
end

topOpt.setFrame( topOpt.findFrame( volumeFractions ) );

% solve natural vibration problem for given volume fraction
vibrations.solve( nEigenForms, topOpt.x );

% saving results in desired folder, separate for each volume fraction
folderName = mainFolderName;

figure;
topOpt.plotCurrentFrame();
savefig( gcf, fullfile( folderName, "self_weight_optimal_topology.fig"));

filename_predictor = fullfile( folderName,"topology_eigen");
fol_fr_text = "Vol. fr. = " + num2str(volumeFractions + " ");
  

multi_forms_results_filename='BeamMultiFormsTopOptConst.mat';

vibrations.plotNaturalForms(folderName+"/self_weight_optimal_topology",1:10,topOpt.x>0.5,"Self weight topology, ",0.02,'k','k');

correlation_matrix = vibrations.ComputeCorrelationMatrix(12,dd_forms(:,1:3))';
vibrations.printCorrelationTable(folderName+"/self_weight_topology_correlation_table","self weigh topology correlation",correlation_matrix);

% If true: always compute topology (ignore existing results)
% If false: compute only if the results file defined in results_filename does not exist
enforce_topology_computation = true;

load_modes_number=3;

if enforce_topology_computation || ~exist(multi_forms_results_filename, "file")

    formTopOpts=[];
    formVibrations=[];

    %topology optimizations for loads proportional to eigenforms
    for k=1:load_modes_number
        analysisHarmonic = ElasticHarmonicVibrations(fe, mesh, k, false, true);
        analysisHarmonic.supports = analysisLinear.supports;
        analysisHarmonic.massClosestNode([l/2 h/2],["ux" "uy"],[10 10])

        % natural vibration analysis object definition
        vibrations = LinearNaturalVibration( analysisLinear.felems, mesh );
        vibrations.supports = analysisLinear.supports;
        vibrations.massClosestNode([l/2 h/2],["ux" "uy"],[10 10])


        topOpt = StressIntensityTopologyOptimizationVol(Rfilter, analysisHarmonic, cutTreshold, penal,  volumeFractions, true);
        topOpt.setConstElems(const_elems);
        figure, hold on;
        [objF, xopt]  = topOpt.solve();
        vibrations.solve( nEigenForms, xopt );
        formVibrations = [formVibrations vibrations];
        formTopOpts = [formTopOpts topOpt];

    end

    vibrations=[];

    adiv=4;
    
    alphas = (1:(adiv-1))/adiv;
    midTopOpts=[];
    frqs=[];

    analysisHarmonic1 = ElasticHarmonicVibrations(fe, mesh, 1, false, true);
    analysisHarmonic1.supports = analysisLinear.supports;
    analysisHarmonic1.massClosestNode([l/2 h/2],["ux" "uy"],[10 10])

    analysisHarmonic2 = ElasticHarmonicVibrations(fe, mesh, 2, false, true);
    analysisHarmonic2.supports = analysisLinear.supports;
    analysisHarmonic2.massClosestNode([l/2 h/2],["ux" "uy"],[10 10])

    analysisHarmonic3 = ElasticHarmonicVibrations(fe, mesh, 3, false, true);
    analysisHarmonic3.supports = analysisLinear.supports;
    analysisHarmonic3.massClosestNode([l/2 h/2],["ux" "uy"],[10 10])

    mixedTopOpts12env=[];
    mixedVibrations12env=[];
    
    % topology optimizations for loads proportional to linear combinations of 1 and 2 eigenforms
    for k=1:adiv-1
        % natural vibration analysis object definition
        vibrations12env = LinearNaturalVibration( analysisLinear.felems, mesh );
        vibrations12env.supports = analysisLinear.supports;

        analysisLinear.selfLoadFactor=-1;
        topOptMulti12 = StressIntensityMultiMaxTopologyOptimization(Rfilter, [analysisHarmonic1 analysisHarmonic2], [ (1-alphas(k)) alphas(k) ], cutTreshold, penal, volumeFractions, true);
        topOptMulti12.setConstElems(const_elems);
        figure, hold on;
        [objF, xopt]  = topOptMulti12.solve();
        vibrations12env.solve( nEigenForms, xopt );
        mixedTopOpts12env = [mixedTopOpts12env topOptMulti12];
        mixedVibrations12env = [mixedVibrations12env vibrations12env];
    end

    mixedTopOpts12av=[];
    mixedVibrations12av=[];
    
    % topology optimizations for loads proportional to linear combinations of 1 and 2 eigenforms
    for k=1:adiv-1
        % natural vibration analysis object definition
        vibrations12av = LinearNaturalVibration( analysisLinear.felems, mesh );
        vibrations12av.supports = analysisLinear.supports;

        analysisLinear.selfLoadFactor=-1;
        topOptMulti12 = StressIntensityMultiAvTopologyOptimization(Rfilter, [analysisHarmonic1 analysisHarmonic2], [ (1-alphas(k)) alphas(k) ], cutTreshold, penal, volumeFractions, true);
        topOptMulti12.setConstElems(const_elems);
        figure, hold on;
        [objF, xopt]  = topOptMulti12.solve();
        vibrations12av.solve( nEigenForms, xopt );
        mixedTopOpts12av = [mixedTopOpts12av topOptMulti12];
        mixedVibrations12av = [mixedVibrations12av vibrations12av];
    end

    mixedTopOpts13env=[];
    mixedVibrations13env=[];
    
    % topology optimizations for loads proportional to linear combinations of 1 and 3 eigenforms
    for k=1:adiv-1
        vibrations13env = LinearNaturalVibration( analysisLinear.felems, mesh );
        vibrations13env.supports = analysisLinear.supports;

        analysisLinear.selfLoadFactor=-1;
        topOptMulti13 = StressIntensityMultiMaxTopologyOptimization(Rfilter, [analysisHarmonic1 analysisHarmonic3], [ (1-alphas(k)) alphas(k) ], cutTreshold, penal, volumeFractions, true);
        topOptMulti13.setConstElems(const_elems);
        figure, hold on;
        [objF, xopt]  = topOptMulti13.solve();
        vibrations13env.solve( nEigenForms, xopt );
        mixedTopOpts13env = [mixedTopOpts13env topOptMulti13];
        mixedVibrations13env = [mixedVibrations13env vibrations13env];
    end

    mixedTopOpts13av=[];
    mixedVibrations13av=[];
    
    % topology optimizations for loads proportional to linear combinations of 1 and 3 eigenforms
    for k=1:adiv-1
        vibrations13av = LinearNaturalVibration( analysisLinear.felems, mesh );
        vibrations13av.supports = analysisLinear.supports;

        analysisLinear.selfLoadFactor=-1;
        topOptMulti13 = StressIntensityMultiAvTopologyOptimization(Rfilter, [analysisHarmonic1 analysisHarmonic3], [ (1-alphas(k)) alphas(k) ], cutTreshold, penal, volumeFractions, true);
        topOptMulti13.setConstElems(const_elems);
        figure, hold on;
        [objF, xopt]  = topOptMulti13.solve();
        vibrations13av.solve( nEigenForms, xopt );
        mixedTopOpts13av = [mixedTopOpts13av topOptMulti13];
        mixedVibrations13av = [mixedVibrations13av vibrations13av];
    end

    % mixedTopOpts23=[];
    % mixedVibrations23=[];
    % 
    % % topology optimizations for loads proportional to linear combinations of 1 and 3 eigenforms
    % for k=1:adiv-1
    %     vibrations23 = LinearNaturalVibration( analysisLinear.felems, mesh );
    %     vibrations23.supports = analysisLinear.supports;
    % 
    %     analysisLinear.selfLoadFactor=-1;
    %     topOptMulti23 = StressIntensityMultiMaxTopologyOptimization(Rfilter, [analysisHarmonic2 analysisHarmonic3], [ (1-alphas(k)) alphas(k) ], cutTreshold, penal, volumeFractions, true);
    %     topOptMulti23.setConstElems(const_elems);
    %     figure, hold on;
    %     [objF, xopt]  = topOptMulti23.solve();
    %     vibrations23.solve( nEigenForms, xopt );
    %     mixedTopOpts23 = [mixedTopOpts23 topOptMulti23];
    %     mixedVibrations23 = [mixedVibrations23 vibrations23];
    % end

    close all;
    save(multi_forms_results_filename,"-v7.3");
    
else
    load(multi_forms_results_filename);
end

results_folder_name = "beam_results_const";

% solve natural vibration problem for design domain
dd_vibrations = LinearNaturalVibration( analysisLinear.felems, mesh );
dd_vibrations.supports = analysisLinear.supports;   % (don't forget this line if needed)
dd_vibrations.solve( nEigenForms, 1 ); 

dd_freq  = dd_vibrations.frequencies;
dd_forms = dd_vibrations.qforms;


% output section

for form_idx=1:3
    
        figure, hold on;
        iters = formTopOpts(form_idx).findFrame( volumeFractions );
        formTopOpts(form_idx).setFrame( iters );
        formTopOpts(form_idx).plotCurrentFrame();
        title("FSD, topology for \Phi_"+num2str(form_idx)+", vol = "+num2str(volumeFractions(1))+", iteration = "+num2str(iters));

        folderName = fullfile(results_folder_name,"topology_for_"+num2str(form_idx)+"_mode");
        if ~isfolder(folderName)
            mkdir(folderName);
        end

        savefig( gcf, fullfile( folderName, "topology_for_"+num2str(form_idx)+"_mode.fig")); 
        formVibrations(form_idx).plotNaturalForms(fullfile(folderName,"topology_for_"+num2str(form_idx)+"_mode"),1:12,formTopOpts(form_idx).x>0.5,", vol =  "+num2str(volumeFractions(1))+", ",0.02,'k','k');
        correlation_matrix = formVibrations(form_idx).ComputeCorrelationMatrix(12,dd_forms(:,1:3))';
        formVibrations(form_idx).printCorrelationTable(fullfile(folderName, "corelations_for_"+num2str(form_idx)+"_mode"),"topology correlation",correlation_matrix);
end

close all;

for k=1:adiv-1
        alpha=alphas(k);
        figure, hold on;
        mixedTopOpts12env(k).setFrame( mixedTopOpts12env(k).findFrame( volumeFractions ));
        mixedTopOpts12env(k).plotCurrentFrame();
        title("FSD, vol="+ num2str(volumeFractions(1))+", P=[ "+num2str(1-alpha)+"*\Phi_1 "+num2str(alpha)+"*\Phi_2]");

        folderName = fullfile(results_folder_name,"topology_for_modes_1_2");
        if ~isfolder(folderName)
            mkdir(folderName);
        end

        file_prefix = "topology_for_modes_1_2_envelope_";

        subFolderName = fullfile(folderName,file_prefix+num2str(k));
        if ~exist(subFolderName, 'dir')
            mkdir(subFolderName);
        end

        savefig( gcf, fullfile( folderName, file_prefix+num2str(k)+".fig")); % % title(["FSD, vol = ..., loading: [0.25*Phi_1  0.75*Phi_3], mode no. ..., freq. = ... Hz"])

        mixedVibrations12env(k).plotNaturalForms(fullfile(subFolderName, file_prefix + "_" + num2str(k)),1:10,mixedTopOpts12env(k).x>0.5," envelope , vol=" + num2str(volumeFractions(1)) + " [ " + num2str(1-alpha)+"*\Phi_1 "+num2str(alpha)+"*\Phi_2 ] ",0.02,'k','k');
        correlation_matrix = mixedVibrations12env(k).ComputeCorrelationMatrix(12,dd_forms(:,1:3))';
        mixedVibrations12env(k).printCorrelationTable(fullfile(subFolderName, file_prefix + "_"+num2str(k)),"topology correlation",correlation_matrix);
end

close all;

for k=1:adiv-1
        alpha=alphas(k);
        figure, hold on;
        mixedTopOpts12av(k).setFrame( mixedTopOpts12av(k).findFrame( volumeFractions ));
        mixedTopOpts12av(k).plotCurrentFrame();
        title("FSD, vol="+ num2str(volumeFractions(1))+", P=[ "+num2str(1-alpha)+"*\Phi_1 "+num2str(alpha)+"*\Phi_2]");

        folderName = fullfile(results_folder_name,"topology_for_modes_1_2");
        if ~exist(folderName, 'dir')
            mkdir(folderName);
        end

        file_prefix = "topology_for_modes_1_2_average_";

        subFolderName = fullfile(folderName,file_prefix+num2str(k));
        if ~isfolder(subFolderName)
            mkdir(subFolderName);
        end

        savefig( gcf, fullfile( folderName, file_prefix+num2str(k)+".fig")); % % title(["FSD, vol = ..., loading: [0.25*Phi_1  0.75*Phi_3], mode no. ..., freq. = ... Hz"])

        mixedVibrations12av(k).plotNaturalForms(fullfile(subFolderName,  file_prefix + "_" + num2str(k)),1:10,mixedTopOpts12av(k).x>0.5," average , vol=" + num2str(volumeFractions(1)) + " [ " + num2str(1-alpha)+"*\Phi_1 "+num2str(alpha)+"*\Phi_2 ] ",0.02,'k','k');
        correlation_matrix = mixedVibrations12av(k).ComputeCorrelationMatrix(12,dd_forms(:,1:3))';
        mixedVibrations12av(k).printCorrelationTable(fullfile(subFolderName,file_prefix + "_"+num2str(k)),"topology correlation",correlation_matrix);
end

close all;

for k=1:adiv-1
        alpha=alphas(k);
        figure, hold on;
        mixedTopOpts13env(k).setFrame( mixedTopOpts13env(k).findFrame( volumeFractions ));
        mixedTopOpts13env(k).plotCurrentFrame();
        title("FSD, vol="+ num2str(volumeFractions(1))+", P=[ "+num2str(1-alpha)+"*\Phi_1 "+num2str(alpha)+"*\Phi_3]");

        folderName = fullfile(results_folder_name,"topology_for_modes_1_3");
        if ~exist(folderName, 'dir')
            mkdir(folderName);
        end

         file_prefix = "topology_for_modes_1_3_envelope_";

        subFolderName = fullfile(folderName,file_prefix+num2str(k));
        if ~isfolder(subFolderName)
            mkdir(subFolderName);
        end

        savefig( gcf, fullfile( folderName, file_prefix + "_" + num2str(k) + ".fig"));

        mixedVibrations13env(k).plotNaturalForms(fullfile(subFolderName,file_prefix+"_"+num2str(k)),1:10,mixedTopOpts13env(k).x>0.5," envelope , vol=" + num2str(volumeFractions(1)) + " [ " + num2str(1-alpha) + "*\Phi_1 "+num2str(alpha)+"*\Phi_3 ] ",0.02,'k','k');
        correlation_matrix = mixedVibrations13env(k).ComputeCorrelationMatrix(12,dd_forms(:,1:3))';
        mixedVibrations13env(k).printCorrelationTable(fullfile(subFolderName,file_prefix+"_"+num2str(k)),"topology correlation",correlation_matrix);
end

close all;

for k=1:adiv-1
        alpha=alphas(k);
        figure, hold on;
        mixedTopOpts13av(k).setFrame( mixedTopOpts13av(k).findFrame( volumeFractions ));
        mixedTopOpts13av(k).plotCurrentFrame();
        title("FSD, vol="+ num2str(volumeFractions(1))+", P=[ "+num2str(1-alpha)+"*\Phi_1 "+num2str(alpha)+"*\Phi_3]");

        folderName = fullfile(results_folder_name,"topology_for_modes_1_3");
        if ~isfolder(folderName)
            mkdir(folderName);
        end

         file_prefix = "topology_for_modes_1_3_average_";

        subFolderName = fullfile(folderName,file_prefix+num2str(k));
        if ~exist(subFolderName, 'dir')
            mkdir(subFolderName);
        end

        savefig( gcf, fullfile( folderName, file_prefix + "_" + num2str(k) + ".fig"));

        mixedVibrations13av(k).plotNaturalForms(fullfile(subFolderName,file_prefix+"_"+num2str(k)),1:10,mixedTopOpts13av(k).x>0.5," average , vol=" + num2str(volumeFractions(1)) + " [ " + num2str(1-alpha) + "*\Phi_1 "+num2str(alpha)+"*\Phi_3 ] ",0.02,'k','k');
        correlation_matrix = mixedVibrations13av(k).ComputeCorrelationMatrix(12,dd_forms(:,1:3))';
        mixedVibrations13av(k).printCorrelationTable(fullfile(subFolderName,file_prefix+"_"+num2str(k)),"topology correlation",correlation_matrix);
end


% for k=1:adiv-1
%         alpha=alphas(k);
%         figure, hold on;
%         mixedTopOpts23(k).setFrame( mixedTopOpts23(k).findFrame( volumeFractions ));
%         mixedTopOpts23(k).plotCurrentFrame();
%         title("FSD P=[ "+num2str(1-alpha)+"*\Phi_2 "+num2str(alpha)+"*\Phi_3]");
% 
%         folderName = fullfile(results_folder_name,"topology_for_modes_2_3");
%         if ~exist(folderName, 'dir')
%             mkdir(folderName);
%         end
% 
%         savefig( gcf, fullfile( folderName, "topology_for_modes_2_3_mix_"+num2str(k)+".fig"));
% 
%         mixedVibrations23(k).plotNaturalForms(folderName+"/topology_for_modes_mix_2_3_form_"+num2str(k),1:10,mixedTopOpts23(k).x>0.5,"topology mix for modes 2 3, "+num2str(1-alpha)+"*\Phi_2 "+num2str(alpha)+"*\Phi_3] ",0.02,'k','k');
%         correlation_matrix = mixedVibrations23(k).ComputeCorrelationMatrix(12,dd_forms(:,1:3))';
%         mixedVibrations23(k).printCorrelationTable(folderName+"/correlations_for_modes_mix_2_3_form_"+num2str(k),"topology correlation",correlation_matrix);
% end


