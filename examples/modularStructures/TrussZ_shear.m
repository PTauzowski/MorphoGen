clear;
close all;


% lobger edge length
l1=1.5;

%shorter edhe length
l2=1;

%height
h=1.5;

%width
w=1;

%thickness
th=0.05;

% FE division along thickness
nth=1;

basename='TrussZ_shear';



%finite element definition
ShapeFn = ShapeFunctionL8;
mesh = Mesh();
mesh.addrectPipe(w,h,l1,th,nth,ShapeFn.localNodes);
fe = SolidElasticElem( ShapeFn, mesh.elems );

%material definition
material = SolidMaterial('mat1');
material.setElasticIzo(1, 0.3);
material.setElasticIzoGrad();
fe.setMaterial(material);

%analysis definition
analysis = LinearElasticityWeighted( fe, mesh, true );

%load forces initialization
N5=0;
N6=0;
Ty5=0;
Ty6=0;
Tz5=0;
Tz6=0;

q=0.1; qy=0.00; lb=1; Tf=0.0;
My=0; Mz=0; Tz=0; Ty=0; Tr=0;

%setting shear force load
%My=q*lb^2/8;     
%Mz=qy*lb^2/8;   
Tz=q*lb/2;      
%Ty=qy*lb/2;     
%Tr=110*q;       


% computation of loads for FE structure
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

% Selectors for selection of fixed nodes
fixedFacesSelector1 = Selector( @(x)( x(:,1)==0 & (x(:,2)<=1.01*th) & (x(:,3)<=1.01*th) ) );
fixedFacesSelector2 = Selector( @(x)( x(:,1)==0 & (w-x(:,2)<=1.01*th) & (x(:,3)<=1.01*th)  ) );
fixedFacesSelector3 = Selector( @(x)( x(:,1)==0 & (w-x(:,2)<=1.01*th) & (h-x(:,3)<=1.01*th ) ) );
fixedFacesSelector4 = Selector( @(x)( x(:,1)==0 & (x(:,2)<=1.01*th) & (h-x(:,3)<=1.01*th ) ) ) ;

% Selectors for selection of const elements (frame which is not subjected
% to elements removal

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

%fixing supported nodes
analysis.fixNodes( fixedFacesSelector1, ["ux" "uy" "uz"] );
analysis.fixNodes( fixedFacesSelector2, ["ux" "uy" "uz"] );
analysis.fixNodes( fixedFacesSelector3, ["ux" "uy" "uz"] );
analysis.fixNodes( fixedFacesSelector4, ["ux" "uy" "uz"] );

%application of nodal load
analysis.loadClosestNode([l1 0, 0],["ux" "uy" "uz"],[N1,Ty1-Tr,Tz1]);
analysis.loadClosestNode([l1 w, 0],["ux" "uy" "uz"],[N2,Ty2,Tz2-Tr]);
analysis.loadClosestNode([l1 w, h],["ux" "uy" "uz"],[N3,Ty3+Tr,Tz3]);
analysis.loadClosestNode([l1 0, h],["ux" "uy" "uz"],[N4,Ty4,Tz4+Tr]);

%transformation of mest to achieve skewnes and differnt height at the
%beginning and end of segment
mesh.transformNodesXY( @(x)( [ x(:,1)+(0.2*(l1-x(:,1))-0.2*x(:,1)).*x(:,2)/l2 x(:,2) x(:,3)+0.3*x(:,1)/l1 ] )  );

%plotting of design space
fe.plotSolid(mesh.nodes);
analysis.plotCurrentLoad();
analysis.plotSupport();
analysis.printProblemInfo();


% Filtering radius
Rfilter = 2*th/nth;

%penalty parameter
penal=3;

% threshold of stress intensity for element removal removal 
cutTreshold = 0.005;

view(45, 45);


tic

% topology otimization object definition
topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penal, 0.45, false );
topOpt.setConstElems(const_elems);

% topology optimization procedure execution
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
