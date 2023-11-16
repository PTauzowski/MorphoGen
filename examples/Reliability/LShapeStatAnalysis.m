clear;
close all;
l=3;
c=0.4;
res = 12;
E=210000;
nu=0.3;
sf = ShapeFunctionL4;
mesh = Mesh();
mesh.addRectMesh2D(0, 0, l, c*l, round(1/c*res), res, sf.pattern);
mesh.addRectMesh2D(0, c*l, c*l, (1-c)*l, res, 1/(1-c)*res, sf.pattern);
fe=PlaneStressElem( sf, mesh.elems );
fe.props.h=1;
material = PlaneStressMaterial('mat1');
material.setElasticIzo(E, nu);
fe.setMaterial( material );

fe.plotSolid(mesh.nodes);
analysis = LinearElasticity( fe, mesh );
fixedEdgeSelectorX = Selector( @(x)( abs(x(:,1))<0.001 ) );
fixedEdgeSelectorY = Selector( @(x)( abs(x(:,2)) )<0.001 );
loadedEdgeSelectorX = Selector( @(x)( abs(x(:,1) - l)<0.001 ) );
loadedEdgeSelectorY = Selector( @(x)( abs(x(:,2) - l)<0.001 ) );

analysis.elementLoadLineIntegral( "global", loadedEdgeSelectorX, ["ux" "uy"], @(x)( x*0 + [ 100 0 ] ));
analysis.elementLoadLineIntegral( "global", loadedEdgeSelectorY, ["ux" "uy"], @(x)( x*0 + [ 0 100 ] ));
analysis.fixNodes( fixedEdgeSelectorX, "ux");
analysis.fixNodes( fixedEdgeSelectorY, "uy" );
analysis.plotCurrentLoad();
analysis.plotSupport();

analysis.printProblemInfo();

figure;
randomVariables={RandomVariable("Normal",100,20) RandomVariable("Normal",100,20)};
g=loadPerformanceFunction(analysis,material,mesh.findClosestNode([c*l c*l]),1,17);
g.loadedEdgeSelectorX=loadedEdgeSelectorX;
g.loadedEdgeSelectorY=loadedEdgeSelectorY;

% N=10000;
% mc= MonteCarlo(randomVariables,g,N);
% [ Pf_mc, p ] = mc.solve();
% Pf_mc
% figure, hold on;
% scatter3(mc.x(p>0,1),mc.x(p>0,2),p(p>0),'MarkerEdgeColor',[0 .8 .8],'Marker','.');
% scatter3(mc.x(p<=0,1),mc.x(p<=0,2),p(p<=0),'filled','MarkerEdgeColor',[0.5 0 .5],'Marker','o');

hmv = HMV(randomVariables,g,3);
form = FORM(randomVariables,g);
Pf_form = form.solve()
[ Pf, mpp, betar ] = hmv.solve();


