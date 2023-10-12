clear;
close all;
res = 8;
l = 1;

Rfilter = 2*l/res;

%Removal intensity threshold
cutTreshold = 0.005;

%penalty factor
penal = 3;

sfL8 = ShapeFunctionL8;
mesh = Mesh();
mesh.addLshape3D( 2*l, 0.8*l, 2*res, sfL8.pattern);
fe=SolidElasticElem( sfL8, mesh.elems );

material = SolidMaterial('mat1');
material.setElasticIzo(1, 0.3);
material.setElasticIzoGrad();
fe.props.h=1;
fe.setMaterial(material)

analysis = LinearElasticityWeighted( fe, mesh, true );
fixedEdgeSelector = Selector( @(x)( x(:,3) - 2*l ) );
analysis.loadClosestNode([ 2*l, 0.4*l, 0.4*l ], ["ux" "uy" "uz"], [0 0 -1.0] );
analysis.fixNodes( fixedEdgeSelector, ["ux" "uy" "uz"] );

analysis.printProblemInfo();
fe.plotSolid(mesh.nodes);
analysis.plotCurrentLoad();
analysis.plotSupport();
 view(45, 45);
 %         view(135, 25);

tic
topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penal, 0.2, true );
[objF, xopt]  = topOpt.solve();
toc

figure;
tic
topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, analysis, penal, 0.2, true);
[objF, xopt]  = topOpt.solve();
toc


