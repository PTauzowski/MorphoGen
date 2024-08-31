clear;
close all;

R=0.1;
th=0.005;
r=R-th;
h=0.2;
alpha=30;

% resCirc=100;
% resTh=max(1,round(th/2/pi/R*resCirc));
% resHeight=max(1,round(h/4/pi/R*resCirc));

resTh=1;
resCirc=round(2*pi*R/th*resTh);
resHeight=round(h/th*resTh);

% Filtering radius
Rfilter = R*3*pi/resCirc;
penal=3;
cutTreshold = 0.005;

ShapeFn = ShapeFunctionL8;
mesh = Mesh();
mesh.addManipulatorHalfSegment3D(r, R, h, alpha*pi/180, resTh, resCirc, resHeight, ShapeFn.localNodes);
fe = SolidElasticElem( ShapeFn, mesh.elems );
fe.plot(mesh.nodes);
upward_facing_nodes=mesh.findUpwardFacingNodes();
downward_facing_nodes=mesh.findDownwardFacingNodes();
% scatter3(mesh.nodes(upward_facing_nodes,1), mesh.nodes(upward_facing_nodes,2), mesh.nodes(upward_facing_nodes,3), ...
%          100, 'r', 'filled');

fe.props.h=1;
material = SolidMaterial('mat1');
material.setElasticIzo(210.0E9, 0.3);
material.setElasticIzoGrad();
fe.setMaterial(material);

analysis = LinearElasticityWeighted( fe, mesh, false );
%problem = LinearElasticity( fe, mesh );
% fixedEdgeSelector = Selector( @(x)( abs(x(:,3)) < 0.001 ) );
% loadedFaceSelector = Selector( @(x)( abs(x(:,3)- Length) < 0.001 ) );
constElemsSelector =  Selector( @(x)( (x(:,3) < 0.05 * h ) ) & (x(:,3) > 0.96 * h ) );

un=false(size(mesh.nodes,1),1);
dn=false(size(mesh.nodes,1),1);
un( upward_facing_nodes ) = true;
dn(downward_facing_nodes ) = true;
loadedFaceSelector = Selector( un );
fixedEdgeSelector = Selector( dn );
upElems = mesh.findElems( loadedFaceSelector, false );
downElems = mesh.findElems( fixedEdgeSelector, false );



%analysis.elementLoadSurfaceIntegral( "global", loadedFaceSelector, ["ux" "uy" "uz"], @(x)( x*0 + [-x(:,2)./sqrt(x(:,1).^2+x(:,2).^2) x(:,1)./sqrt(x(:,1).^2+x(:,2).^2) -x(:,2)./x(:,2)] )); 
analysis.elementLoadSurfaceIntegral( "global", loadedFaceSelector, ["ux" "uy" "uz"], @(x)( x*0 + [0 0 -2.0E8] ));
analysis.fixNodes( fixedEdgeSelector, ["ux" "uy" "uz"] );
%analysis.fixClosestNode( [0 0 0], ["ux" "uy" "uz"], [0 0 0]);
const_elems = [upElems downElems];

%mesh.transformNodesXY( @(x)( [ x(:,1) x(:,2) x(:,3)-0.3*x(:,1).*x(:,3)/Length ] )  );

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

analysisSecondOrder = SecondOrderElasticityWeighted( fe, mesh, 0.0, false );
analysisSecondOrder.Pnodal=stability.Pnodal;
analysisSecondOrder.Pfem=stability.Pfem;
analysisSecondOrder.supports=stability.supports;

analysisWithBuckling = SecondOrderElasticityWeighted( fe, mesh, 0.95, false );
analysisWithBuckling.Pnodal=stability.Pnodal;
analysisWithBuckling.Pfem=stability.Pfem;
analysisWithBuckling.supports=stability.supports;

% figure;
% tic
% topOptLinear = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penal, 0.2, false );
% [objF, xopt]  = topOptLinear.solve();
% toc

figure;
tic
topOptSecondOrder = StressIntensityTopologyOptimizationBuckling( Rfilter, analysisSecondOrder, cutTreshold, penal, 0.4, true );
topOptSecondOrder.setConstElems(const_elems);
[objF, xopt]  = topOptSecondOrder.solve();
toc

figure;
tic
topOptBuckling = StressIntensityTopologyOptimizationBuckling( Rfilter, analysisWithBuckling, cutTreshold, penal, 0.4, true );
topOptBuckling.setConstElems(const_elems);
[objF, xopt]  = topOptBuckling.solve();
toc


