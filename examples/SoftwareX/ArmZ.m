clear;
close all;

width=200;
length=250;
height=100;
r=40;
R=70;
xs1=[width/2 0];
xs2=[width/2 length];
alpha=deg2rad(30);
m=tan(alpha);
c1=xs1(2)-m*xs1(1);
c2=xs2(2)+m*xs2(1);
resC=2;

cp=0.001;
Ps=cp*4;
Ph=cp*1;

sfL8= ShapeFunctionL8();
mesh = snakeSegmentModel(width, length, height, r, R, sfL8.localNodes,resC);
halfSymetricalSelector = Selector( @(x)(  (width/2 - x(:,1)  ) > 0.01 ) );
mesh.removeNodes( halfSymetricalSelector );

fe = SolidElasticElem( sfL8, mesh.elems );
maxMesh = max(mesh.nodes);
minMesh = min(mesh.nodes);
l=max(maxMesh-minMesh);

x=ones(1,size(mesh.elems,1));

fe.props.h=1;
material = SolidMaterial('mat1');
material.setElasticIzo(1, 0.3);
material.setElasticIzoGrad();
fe.setMaterial(material);

%problem = LinearElasticity( fe, mesh );
problem = LinearElasticityWeighted( fe, mesh, false );
fixedFaceSelector = Selector( @(x)( abs(x(:,3))<0.001 ) );
fixedSymetrySelector = Selector( @(x)( abs(x(:,1) - width/2) < 0.001 ) );

loadedFaceSelector = Selector( @(x)( abs( (x(:,1) - xs2(1)).^2  + (x(:,2) - xs2(2)).^2 - r^2) < 0.0001 ) );
  holeFaceSelector = Selector( @(x)( abs( (x(:,1) - xs1(1)).^2  + (x(:,2) - xs1(2)).^2 - r^2) < 0.0001 ) );

sideFaceSelector1 = Selector( @(x)( abs(m*x(:,1)+c1-x(:,2)) < 0.01) );
sideFaceSelector2 = Selector( @(x)( abs(-m*x(:,1)+c2-x(:,2)) < 0.01) );

problem.fixNodes( fixedFaceSelector, ["uz"] );
problem.fixNodes( holeFaceSelector, ["ux" "uy"] );
problem.fixNodes( fixedSymetrySelector, ["ux"] );
problem.elementLoadSurfaceIntegral( "global", loadedFaceSelector, ["ux" "uy" "uz"], @(x)( x.*0 + [0 Ph 0] ));
problem.elementLoadSurfaceIntegral( "global", sideFaceSelector1, ["ux" "uy" "uz"], @(x)( x.*0 + [-Ps*sin(alpha)   Ps*cos(alpha) 0] ));
problem.elementLoadSurfaceIntegral( "global", sideFaceSelector2, ["ux" "uy" "uz"], @(x)( x.*0 + [-Ps*sin(alpha)  -Ps*cos(alpha) 0] ));


fe.plotSolid(mesh.nodes);
problem.printProblemInfo();
problem.plotCurrentLoad();
problem.plotSupport();
view(45, 45);

% Filtering radius
Rfilter = 2*l/(resC*32);
penal=3;
cutTreshold = 0.005;


tic
topOpt = StressIntensityTopologyOptimizationVol( Rfilter, problem, cutTreshold, penal, 0.3, false );
[objF, xopt]  = topOpt.solve();
toc
save('ArmZ.mat');

% tic
% topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, problem, penal, 0.3, false);
% [objF, xopt]  = topOpt.solve();
