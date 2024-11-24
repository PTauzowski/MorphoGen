
clear;
close all;

Ecl  = 25.0E9; 
Ebm  = 25.0E9;
Egr  = 80.0E9;
nugr = 0.3;

b=0.25;
hbm=0.8;
hcl=0.7;
lspan=5;
hfloor=3;
nspan=3;
nfloor=3;


lgr=30;
hgr=15;
resgr=5;

%model = FrameOnElasticGroundModel( Ecl, hcl, Ebm, hbm, b, lspan, hfloor, nspan, nfloor, Egr, nugr, lgr, hgr, resgr, ShapeFunctionL4 );


xp=[(lgr-nspan*lspan)/2 hgr];

shapeFn2D=ShapeFunctionL16();
       
mesh=Mesh();
frameElem = Frame2D( mesh.addHframe(nspan,lspan,nfloor,hfloor,xp));
planeElem = PlaneStressElem( shapeFn2D, mesh.addRectMeshArray2D(0,0,[xp(1) repelem(lspan,1,nspan) xp(1)], xp(2), resgr, shapeFn2D.pattern) );

model = FEModel( { frameElem planeElem }, mesh );
te=model.getFEClassesNumber();
model.draw();
edgeX0 = selectX(0.0);
model.fixDOF( edgeX0, ["ux" "uy"] , [0 0]);
%model.plotSelectedNodes(edgeX0,'m');
%model.plotNodes(".",'b');

% Example Input
A = {1:3; [4 5]};
B = {2:4; [5 6]};

% Create C where each row is the union of corresponding rows from A and B
C = cellfun(@union, A, B, 'UniformOutput', false);

% Display result
disp(C);
% Output:
% { [1 2 3 4]      % Union of A{1} and B{1}
%   [4 5 6] }      % Union of A{2} and B{2}

