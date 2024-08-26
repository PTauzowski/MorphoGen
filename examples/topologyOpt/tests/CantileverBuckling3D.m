clear;
close all;
res = 20;
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
fixedFaceSelector = Selector( @(x)( abs( x(:,1)) < 0.001 ) );
symmetrySelector = Selector( @(x)( abs(x(:,2) - l) < 0.001 ) );
analysis.loadClosestNode( [2*l l 0], "uz", -1 );
analysis.fixNodes( fixedFaceSelector, ["ux" "uy" "uz"] );
analysis.fixNodes( symmetrySelector, ["uy"] );
%problem.fixClosestNode([0 0 0],["ux" "uy" "uz"],[0 0 0])

analysis.printProblemInfo();
fe.plot(mesh.nodes);
analysis.plotCurrentLoad();
analysis.plotSupport();

nEigenForms=10;
stability = LinearStability( analysis.felems, mesh);
stability.Pnodal = analysis.Pnodal;
stability.Pfem = analysis.Pfem;
stability.supports = analysis.supports;
stability.solve( nEigenForms);
lambdas = diag(stability.lambdas)
for k=1:min(10,nEigenForms)
    subplot(5, 2, k);
    stability.setForm(k);
    fe.plotSolidDeformed(mesh.nodes,stability.qnodal,0.2);
    axis on, xlabel('x-axis'), ylabel('y-axis'), view(3)
    lambda_str = sprintf('%.4g', lambdas(k));
    title(['Form:' num2str(k), ' \lambda=' lambda_str]);
end

%analysis = SecondOrderElasticityWeighted( fe, mesh, lambdas(nEigenForms), false );
analysis = LinearElasticityWeighted( fe, mesh, false );
analysis.Pnodal=stability.Pnodal;
analysis.Pfem=stability.Pfem;
analysis.supports=stability.supports;

figure;
tic
topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penal, 0.2, true );
[objF, xopt]  = topOpt.solve();
toc

% figure;
% tic
% topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, analysis, penal, 0.2, true);
% [objF, xopt]  = topOpt.solve();
% toc


