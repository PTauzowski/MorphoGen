clear;
close all;
res = 32;
l = 3;

% Filtering radius
Rfilter = 2*l/res;
penal=3;
cutTreshold = 0.005;


ShapeFn = ShapeFunctionL8;
mesh = Mesh();
mesh.addRectMesh3D( 0, 0, 0, 2*l, l, l, 2*res, res, res, ShapeFn.localNodes);
fe = SolidElasticElem( ShapeFn, mesh.elems );
mX = max(mesh.nodes(:,1));

fe.props.h=1;
material = SolidMaterial('mat1');
material.setElasticIzo(1, 0.3);
material.setElasticIzoGrad();
fe.setMaterial(material)

analysis = LinearElasticityWeighted( fe, mesh, true );
fixedFaceSelector = Selector( @(x)( x(:,1) ) );
loadedFaceSelector = Selector( @(x)( x(:,1) - 2*l ) );
symmetrySelector = Selector( @(x)( x(:,2) - l ) );
circleSelector = Selector( @(x)( ((x(:,1) - 1.5).^2 + (x(:,2) - 1.5).^2 )-1.5 ), 0.2 );
analysis.loadClosestNode( [2*l l 0], "uz", -1 );
analysis.fixNodes( fixedFaceSelector, ["ux" "uy" "uz"] );
analysis.fixNodes( symmetrySelector, ["uy"] );
%problem.fixClosestNode([0 0 0],["ux" "uy" "uz"],[0 0 0])

analysis.printProblemInfo();
fe.plotSolid(mesh.nodes);
analysis.plotCurrentLoad();
analysis.plotSupport();

tic
topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penal, 0.2, true );
[objF, xopt]  = topOpt.solve();
toc

figure;
tic
topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, analysis, penal, 0.2, true);
[objF, xopt]  = topOpt.solve();
toc


