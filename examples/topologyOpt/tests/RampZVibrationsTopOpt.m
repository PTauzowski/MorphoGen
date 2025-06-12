clear;
clear classes
clear all
rehash
close all;

% Cantilever topology optimization elastic task.

% Resolution of shortest (vertical) edge
res = 50;

% height of the cantilever
h = 1;

% Aspect ratio length/height
aspect=12;

% length of the cantilever
l=aspect*h;

% Filtering radius
Rfilter = 4*h/res;

%Removal intensity threshold
cutTreshold = 0.005;

%penalty factor
penal = 3;

% Type of shape function to be used (here: four node Langrange)
sfL4 = ShapeFunctionL4;

% Creating FE mesh object
mesh = Mesh();

% Generating rectangular mesh ( aspect*h x h )
mesh.addRectMesh2D(0, 0, l, h, aspect*res, res, sfL4.pattern);

% Create plane stress finite element object
fe=PlaneStressElem( sfL4, mesh.elems );

% Create isotropic material object
material = PlaneStressMaterial('mat1');
material.setElasticIzo(205E9, 0.3);
material.setMassIzoMatrix(7850);
% Assigning material to finite element
fe.setMaterial( material );

% Ke = fe.computeStifnessMatrix(mesh.nodes,1,1);
% Me = fe.computeMassMatrix(mesh.nodes,1,1);
% 
%  [q, l] = eigs( Ke(:,:,1),Me(:,:,1), 8, 'smallestabs');

% Creating linear elastic finite element analysis object with weighted matrix feature, weighted by element density.
analysisLinear = LinearElasticityWeighted( fe, mesh, false );

% Creating node selector object to select fixed edge (left)
loadEdgeSelector = Selector( @(x)( abs(x(:,2) -h ) < 0.0005 ) );
fixedEdgeSelector = Selector( @(x)( abs(x(:,1)) < 0.0005 ) );
fixedEdgeSelector2 = Selector( @(x)( abs(x(:,1)-l) < 0.0005 ) );

% Fixing structure according to above defined node selector object
analysisLinear.fixNodes( fixedEdgeSelector, ["ux" "uy"] );
analysisLinear.fixNodes( fixedEdgeSelector2, ["ux" "uy"] );

% Creating load vector with one node loaded at the middle of right edge
P=-2.0E9; %100;
%P=-1.5E8; %150;

hp=h;
analysisLinear.loadClosestNode([ l/2, hp ], ["ux" "uy"], [0 P] );
%analysisLinear.elementLoadLineIntegral( "global",loadEdgeSelector, ["ux" "uy"], @(x)( x*0 + [0 P/l] ));


const_rows=3;
ncel=round(const_rows*res*l);
const_elems=[1:ncel size(mesh.elems,1):-1:size(mesh.elems,1)-ncel ];

% topOpt=FreeVibrationsTopologyOpt(analysisLinear,mesh,'FreeVib');
% topOpt.solve()

nEigenForms=10;
vibrations = LinearNaturalVibration( analysisLinear.felems, mesh );
vibrations.Pnodal = analysisLinear.Pnodal;
vibrations.Pfem = analysisLinear.Pfem;
vibrations.supports = analysisLinear.supports;
vibrations.solve( nEigenForms);

harmonicVivrations = ElasticHarmonicVibrations(analysisLinear.felems, mesh, 1, true, true);
harmonicVivrations.Pnodal = analysisLinear.Pnodal;
harmonicVivrations.Pfem = analysisLinear.Pfem;
harmonicVivrations.supports = analysisLinear.supports;

nforms=3;

basename='const_load';

scales = [0.05 0.05 0.05; 0.05 0.05 0.05; 0.05 0.05 0.05];
for k=1:nforms
    harmonicVivrations.count=0;
    harmonicVivrations.mode=k;
    topOptHarmonicVivrations = StressIntensityTopologyOptimizationDynamicBuckling( Rfilter, harmonicVivrations, cutTreshold, penal, 0.45, false );
    topOptHarmonicVivrations.setConstElems(const_elems);
    [objF, xopt]  = topOptHarmonicVivrations.solve();
    save([basename '_' num2str(k) '.mat'], '-v7.3');
    topOptHarmonicVivrations.plot_frequencies(basename, nforms)
    topOptHarmonicVivrations.plot_forms(basename, nforms, topOptHarmonicVivrations.iteration-1, scales(k,:) );
    topOptHarmonicVivrations.plot_correlation_curve( basename, nforms)
    topOptHarmonicVivrations.plot_uncorrelated_frames( basename )
    topOptHarmonicVivrations.resetAnalysis();
    harmonicVivrations.modes=[];
end

% harmonicVivrations.isLoadConst=false;
% 
% basename='Var_load';
% for k=1:nforms
%     harmonicVivrations.count=0;
%     harmonicVivrations.mode=k;
%     topOptHarmonicVivrations = StressIntensityTopologyOptimizationDynamicBuckling( Rfilter, harmonicVivrations, cutTreshold, penal, 0.45, false );
%     topOptHarmonicVivrations.setConstElems(const_elems);
%     [objF, xopt]  = topOptHarmonicVivrations.solve();
%     topOptHarmonicVivrations.plot_frequencies(basename, nforms)
%     topOptHarmonicVivrations.plot_forms(basename, nforms, topOptHarmonicVivrations.iteration-1)
%     topOptHarmonicVivrations.plot_correlation_curve( basename, nforms)
%     topOptHarmonicVivrations.plot_uncorrelated_frames( basename)
%     topOptHarmonicVivrations.resetAnalysis();
%     harmonicVivrations.modes=[];
% end





%save('Cantilever2DBucklingDown80.mat');

