clear;
close all;

R=5;
Th=0.2;
Length=30;

resLen=100;
resCirc = ceil(resLen/Length*2*pi*R);
resTh = ceil(resLen/Length*Th);


% Filtering radius
Rfilter = R*3*pi/resCirc;
penal=3;
cutTreshold = 0.05;

ShapeFn = ShapeFunctionL8;
mesh = Mesh();
mesh.addRectMesh3D( R-Th, 0, 0, Th, 2*pi, Length, resTh, resCirc, resLen, ShapeFn.localNodes);
mesh.transformToCylindrical3D( [0 0] );
fe = SolidElasticElem( ShapeFn, mesh.elems );

fe.props.h=1;
material = SolidMaterial('mat1');
material.setElasticIzo(1, 0.3);
material.setElasticIzoGrad();
fe.setMaterial(material);

analysis = LinearElasticityWeighted( fe, mesh, false );
%problem = LinearElasticity( fe, mesh );
fixedEdgeSelector = Selector( @(x)( abs(x(:,3)) < 0.001 ) );
loadedFaceSelector = Selector( @(x)( abs(x(:,3)- Length) < 0.001 ) );
constElemsSelector = @(x)( (x(:,3) < 0.05 * Length ) ) & (x(:,3) > 0.96 * Length );

analysis.elementLoadSurfaceIntegral( "global", loadedFaceSelector, ["ux" "uy" "uz"], @(x)( x*0 + [-x(:,2)./sqrt(x(:,1).^2+x(:,2).^2) x(:,1)./sqrt(x(:,1).^2+x(:,2).^2) -x(:,2)./x(:,2)] ));
analysis.fixNodes( fixedEdgeSelector, ["ux" "uy" "uz"] );
analysis.fixClosestNode( [0 0 0], ["ux" "uy" "uz"], [0 0 0]);
const_elems = analysis.selectElems( constElemsSelector );

mesh.transformNodesXY( @(x)( [ x(:,1) x(:,2) x(:,3)-0.3*x(:,1).*x(:,3)/Length ] )  );

fe.plotSolid(mesh.nodes);
%problem.plotNodes();
analysis.plotCurrentLoad();
analysis.plotSupport();

view(45, 45);

analysis.printProblemInfo();

tic
topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penal, 0.4, true );
[objF, xopt]  = topOpt.solve();
toc

figure;
tic
topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, analysis, penal, 0.4, true);
[objF, xopt]  = topOpt.solve();
toc




