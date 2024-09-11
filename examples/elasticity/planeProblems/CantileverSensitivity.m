clear;
close all;
res = 10;
l = 3;
tic
sfL4 = ShapeFunctionL16;
mesh = Mesh();
elems = mesh.addRectMesh2D(0, 0, 2*l, l, 2*res, res, sfL4.pattern);
fixedEdgeSelector = Selector( @(x)( abs(x(:,1))<0.001 ) );
loadedEdgeSelector = Selector( @(x)( abs(x(:,1) - 2*l)<0.001 ) );

fe=PlaneStressElem( sfL4, elems );
material = PlaneStressMaterial('mat1');
material.setElasticIzo(1, 0.3);
material.setElasticIzoGrad();
fe.setMaterial( material );
fe.props.h=1;
fe.plot(mesh.nodes);

analysis = LinearElasticityWithMatSensitivity( fe, elems );
analysis.elementLoadLineIntegral( "global",loadedEdgeSelector, ["ux" "uy"], @(x)( x*0 + [0 150] ));
analysis.fixNodes( fixedEdgeSelector, ["ux" "uy"] );
analysis.plotCurrentLoad();
analysis.plotSupport();
analysis.printProblemInfo();

[q, dq] = analysis.solveWithGradient();

eps=0.0001;
material2 = PlaneStressMaterial('mat2');
material2.setElasticIzo(1+eps, 0.3);
fe.setMaterial( material2 );
q2 = analysis.solve();

material3 = PlaneStressMaterial('mat3');
material3.setElasticIzo(1, 0.3+eps);
fe.setMaterial( material3 );
q3  = analysis.solve();

%problem.plotMaps(q,["ux" "uy" "sxx" "syy" "sHM"],0.1);
fe.plotWired(mesh.nodes,q,0.1);

fe.setMaterial( material );

tic

toc

dqA=[ analysis.toFEMVector(q2)-analysis.toFEMVector(q) analysis.toFEMVector(q3)-analysis.toFEMVector(q)]./eps;

disp('');
disp('Exac sensitivity');
max(dq)

disp('Finite difference method');
max(dqA)


fnG=@(x)( sin(x) );
gradG=@(x)(cos(x)  );

%[an,num,error] = sensitivityCheck(fnG,gradG,0.1);

fnG=@(x)( [sin(x(1))*cos(x(2)) cos(x(1))*sin(x(2))] );
gradG=@(x)( [ cos(x(1))*cos(x(2)) -sin(x(1))*sin(x(2)); -sin(x(1))*sin(x(2)) cos(x(1))*cos(x(2))]  );

%[an,num,error] = sensitivityCheck(fnG,gradG,[0.1 0.1]);

