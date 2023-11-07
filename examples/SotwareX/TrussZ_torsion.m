clear;
close all;

l1=1.5;
l2=1;
h=1.5;
w=1;
th=0.05;
nth=1;

basename='TrussZ_torsion';


N1=0.0;
N2=0.0;
N3=0.0;
N4=0.0;
Ty1=0;
Ty2=0;
Ty3=0;
Ty4=0;
Tz1=0;
Tz2=0;
Tz3=0;
Tz4=0;

N5=0;
N6=0;
Ty5=0;
Ty6=0;
Tz5=0;
Tz6=0;

q=0.1; qy=0.00; lb=1; Tf=0.0;

My=q*lb^2/8;     My=0;
Mz=qy*lb^2/8;    Mz=0;
Tz=q*lb/2;       Tz=0;
Ty=qy*lb/2;      Ty=0;
Tr=110*q;        %Tr=0;

N1=My/h/th/th;
N2=My/h/th/th;
N3=-My/h/th/th;
N4=-My/h/th/th;

N1=N1+Mz/h/th/th;
N2=N2+Mz/h/th/th;
N3=N3-Mz/h/th/th;
N4=N4-Mz/h/th/th;

Tz1=Tz/2/th/th;
Tz2=Tz/2/th/th;
Tz3=Tz/2/th/th;
Tz4=Tz/2/th/th;

Ty1=Ty/2/th/th;
Ty2=Ty/2/th/th;
Ty3=Ty/2/th/th;
Ty4=Ty/2/th/th;

ShapeFn = ShapeFunctionL8;
mesh = Mesh();
mesh.addrectPipe(w,h,l1,th,nth,ShapeFn.localNodes);
fe = SolidElasticElem( ShapeFn, mesh.elems );

fe.props.h=1;
material = SolidMaterial('mat1');
material.setElasticIzo(1, 0.3);
material.setElasticIzoGrad();
fe.setMaterial(material);
analysis = LinearElasticityWeighted( fe, mesh, true );

fixedFacesSelector1 = Selector( @(x)( x(:,1)==0 & (x(:,2)<=1.01*th) & (x(:,3)<=1.01*th) ) );
fixedFacesSelector2 = Selector( @(x)( x(:,1)==0 & (w-x(:,2)<=1.01*th) & (x(:,3)<=1.01*th)  ) );
fixedFacesSelector3 = Selector( @(x)( x(:,1)==0 & (w-x(:,2)<=1.01*th) & (h-x(:,3)<=1.01*th ) ) );
fixedFacesSelector4 = Selector( @(x)( x(:,1)==0 & (x(:,2)<=1.01*th) & (h-x(:,3)<=1.01*th ) ) ) ;

constElemsSelectorGate1 = Selector( @(x)( (x(:,1)<=1.1*th) ) ) ;
constElemsSelectorGate2 = Selector( @(x)( (l1-x(:,1)<=1.1*th) ) );

constBarSelector1 = Selector( @(x)( (x(:,1)>0.95*th) & ((l1-x(:,1))>0.95*th) &  (x(:,2)<=1.01*th) & (x(:,3)<1.01*th ) ) );
constBarSelector2 = Selector( @(x)( (x(:,1)>0.95*th) & ((l1-x(:,1))>0.95*th) &  ((w-x(:,2))<=1.01*th) & (x(:,3)<1.01*th)  ) );
constBarSelector3 = Selector( @(x)( (x(:,1)>0.95*th) & ((l1-x(:,1))>0.95*th) &  ((w-x(:,2))<=1.01*th) & ((h-x(:,3))<1.01*th ) ) );
constBarSelector4 = Selector( @(x)( (x(:,1)>0.95*th) & ((l1-x(:,1))>0.95*th) &  (x(:,2)<=1.01*th) & ((h-x(:,3))<1.01*th ) ) ) ;

floorFacesSelector = Selector( @(x)( abs(x(:,3)-th)<0.001 ) );

const_elems = [ mesh.findElems( constElemsSelectorGate1 ); ...
                mesh.findElems( constElemsSelectorGate2 ); ...
                mesh.findElems( constBarSelector1 ); ...
                mesh.findElems( constBarSelector2 );...
                mesh.findElems( constBarSelector3 );...
                mesh.findElems( constBarSelector4 )];

analysis.fixNodes( fixedFacesSelector1, ["ux" "uy" "uz"] );
analysis.fixNodes( fixedFacesSelector2, ["ux" "uy" "uz"] );
analysis.fixNodes( fixedFacesSelector3, ["ux" "uy" "uz"] );
analysis.fixNodes( fixedFacesSelector4, ["ux" "uy" "uz"] );


loadedFacesSelector1 = Selector( @(x)( abs(x(:,1)-l1)<0.001 &  (x(:,2)<=1.1*th) & (x(:,3)<=1.1*th) ) );
loadedFacesSelector2 = Selector( @(x)( abs(x(:,1)-l1)<0.001 &  (w-x(:,2)<=1.1*th) & (x(:,3)<=1.1*th)  ) );
loadedFacesSelector3 = Selector( @(x)( abs(x(:,1)-l1)<0.001 &  (w-x(:,2)<=1.1*th) & (h-x(:,3)<=1.1*th ) ) );
loadedFacesSelector4 = Selector( @(x)( abs(x(:,1)-l1)<0.001 &  (x(:,2)<=1.1*th) & (h-x(:,3)<=1.1*th )  )) ;
loadedFacesSelector5 = Selector( @(x)( abs(x(:,1)-l1)<0.001 &  (x(:,2)<=1.1*th) & ( abs(x(:,3)-h/2)<=1.1*th) ) );
loadedFacesSelector6 = Selector( @(x)( abs(x(:,1)-l1)<0.001 &  (w-x(:,2)<=1.1*th) & (abs(x(:,3)-h/2)<=1.1*th)  ) );

analysis.loadClosestNode([l1 0, 0],["ux" "uy" "uz"],[N1,Ty1-Tr,Tz1]);
analysis.loadClosestNode([l1 w, 0],["ux" "uy" "uz"],[N2,Ty2,Tz2-Tr]);
analysis.loadClosestNode([l1 w, h],["ux" "uy" "uz"],[N3,Ty3+Tr,Tz3]);
analysis.loadClosestNode([l1 0, h],["ux" "uy" "uz"],[N4,Ty4,Tz4+Tr]);

mesh.transformNodesXY( @(x)( [ x(:,1)+(0.2*(l1-x(:,1))-0.2*x(:,1)).*x(:,2)/l2 x(:,2) x(:,3)+0.3*x(:,1)/l1 ] )  );
fe.plotSolid(mesh.nodes);

analysis.plotCurrentLoad();
analysis.plotSupport();
analysis.printProblemInfo();


% Filtering radius
Rfilter = 2*th/nth;
penal=3;
cutTreshold = 0.005;

view(45, 45);

% xall=ones(size(mesh.elems,1));
% tic
% q = analysis.solveWeighted(xall);
% analysis.computeElementResults(q,xall);
% toc
% 
% analysis.plotMaps([ "sxx" "szz" "sHM"],0.1);
% fe.plotWired(mesh.nodes,analysis.qnodal,0.1);

tic
topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penal, 0.45, false );
topOpt.setConstElems(const_elems);
[objF, xopt]  = topOpt.solve();
toc;
save([basename '.mat']);

% 
% figure;
% tic
% topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, analysis, penal, 0.4, false);
% topOpt.setConstElems(const_elems);
% [objF, xopt]  = topOpt.solve();


toc
