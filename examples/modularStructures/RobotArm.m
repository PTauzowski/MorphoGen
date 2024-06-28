clear;
close all;
R=1;
r=0.8;
l=[2 2 2 2];
angles=[30 30 30 -30];
angle=20;
res=3;

Rfilter = 3*(R-r)/res;
%Removal intensity threshold
cutTreshold = 0.005;

%penalty factor
penal = 3;


%finite element definition
ShapeFn = ShapeFunctionL8;
mesh = Mesh();
mesh.addRobotArm(r,R,res,angle,angles,l,ShapeFn.localNodes);
fe = SolidElasticElem( ShapeFn, mesh.elems );
fe.plotSolid(mesh.nodes);

%material definition
material = SolidMaterial('mat1');
material.setElasticIzo(1, 0.3);
material.setElasticIzoGrad();
fe.setMaterial(material);

%analysis definition
analysis = LinearElasticityWeighted( fe, mesh, true );

fixedEdgeSelector = Selector( @(x)( abs(x(:,3)) < 0.0001 ) );
loadedEdgeSelector = Selector( @(x)( (abs( x(:,1) - aspect*l) < 1.0E-4) & ( abs(x(:,2)-l)<l/15 )  ) );

[v,i]=max(mesh.nodes(:,3));

analysis.loadClosestNode(mesh.nodes(i,:), ["ux" "uy" "uz"], [0 0 -1] );

% analysis.plotCurrentLoad();
% analysis.plotSupport();

analysis.printProblemInfo();


% tic
% topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penal, 0.4, false );
% [objF, xopt]  = topOpt.solve();
% toc
% 
% figure;
% tic
% topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, analysis, penal, 0.4, false);
% [objF, xopt]  = topOpt.solve();
% toc
