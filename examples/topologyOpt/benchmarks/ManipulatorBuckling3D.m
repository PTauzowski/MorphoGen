clear;
close all;

R=5;
Th=0.2;
Length=10;

resLen=50;
resCirc = ceil(resLen/Length*2*pi*R);
resTh = ceil(resLen/Length*Th);


% Filtering radius
Rfilter = R*3*pi/resCirc;
penal=3;
cutTreshold = 0.005;

ShapeFn = ShapeFunctionL8;
mesh = Mesh();
mesh.addRectMesh3D( R-Th, 0, 0, Th, 2*pi, Length, resTh, resCirc, resLen, ShapeFn.localNodes);
mesh.transformToCylindrical3D( [0 0] );
fe = SolidElasticElem( ShapeFn, mesh.elems );

fe.props.h=1;
material = SolidMaterial('mat1');
material.setElasticIzo(210.0E9, 0.3);
material.setElasticIzoGrad();
fe.setMaterial(material);

analysis = LinearElasticityWeighted( fe, mesh, false );
%problem = LinearElasticity( fe, mesh );
fixedEdgeSelector = Selector( @(x)( abs(x(:,3)) < 0.001 ) );
loadedFaceSelector = Selector( @(x)( abs(x(:,3)- Length) < 0.001 ) );
constElemsSelector =  Selector( @(x)( (x(:,3) < 0.05 * Length ) ) & (x(:,3) > 0.96 * Length ) );

%analysis.elementLoadSurfaceIntegral( "global", loadedFaceSelector, ["ux" "uy" "uz"], @(x)( x*0 + [-x(:,2)./sqrt(x(:,1).^2+x(:,2).^2) x(:,1)./sqrt(x(:,1).^2+x(:,2).^2) -x(:,2)./x(:,2)] )); 
analysis.elementLoadSurfaceIntegral( "global", loadedFaceSelector, ["ux" "uy" "uz"], @(x)( x*0 + [0 0 -1.0E9] ));
analysis.fixNodes( fixedEdgeSelector, ["ux" "uy" "uz"] );
%analysis.fixClosestNode( [0 0 0], ["ux" "uy" "uz"], [0 0 0]);
const_elems = analysis.selectElems( constElemsSelector );

mesh.transformNodesXY( @(x)( [ x(:,1) x(:,2) x(:,3)-0.3*x(:,1).*x(:,3)/Length ] )  );

fe.plot(mesh.nodes);
%problem.plotNodes();
analysis.plotCurrentLoad();
analysis.plotSupport();

view(45, 45);

analysis.printProblemInfo();

nEigenForms=10;
stability = LinearStability( analysis.felems, mesh);
stability.Pnodal = analysis.Pnodal;
stability.Pfem = analysis.Pfem;
stability.supports = analysis.supports;
stability.solve( nEigenForms);
lambdas = diag(stability.lambdas)
for k=1:min(10,nEigenForms)
    %figure;
    subplot(5, 2, k);
    stability.setForm(k);
    fe.plotSolidDeformed(mesh.nodes,stability.qnodal,0.2);
    axis on, xlabel('x-axis'), ylabel('y-axis'), view(3)
    lambda_str = sprintf('%.4g', lambdas(k));
    title(['Form:' num2str(k), ' \lambda=' lambda_str]);
end

analysisWithBuckling = SecondOrderElasticityWeighted( fe, mesh, 0.90, false );
analysisWithBuckling.Pnodal=stability.Pnodal;
analysisWithBuckling.Pfem=stability.Pfem;
analysisWithBuckling.supports=stability.supports;

analysisWithBuckling = SecondOrderElasticityWeighted( fe, mesh, 0.90, false );
analysisWithBuckling.Pnodal=stability.Pnodal;
analysisWithBuckling.Pfem=stability.Pfem;
analysisWithBuckling.supports=stability.supports;

figure;
tic
topOptLinear = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penal, 0.2, false );
[objF, xopt]  = topOptLinear.solve();
toc

% figure;
% tic
% topOptSecondOrder = StressIntensityTopologyOptimizationBuckling( Rfilter, analysisSecondOrder, cutTreshold, penal, 0.38, true );
% [objF, xopt]  = topOptSecondOrder.solve();
% toc
% 
% figure;
% tic
% topOptBuckling = StressIntensityTopologyOptimizationBuckling( Rfilter, analysisWithBuckling, cutTreshold, penal, 0.38, true );
% [objF, xopt]  = topOptBuckling.solve();
% toc


