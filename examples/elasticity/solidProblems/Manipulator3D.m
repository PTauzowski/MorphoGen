 
clear;
close all;

R=0.5;
r=0.4;
res=20;
alpha=30;
betas=[ 23 15 83 42 56 ];
ShapeFn = ShapeFunctionL8;
segmentLength=0.6;
model = ManipulatorModel3D(segmentLength,R,r,res, alpha, betas, ShapeFn);

model.fe.plot(model.mesh.nodes);
np=100;

% betas=random('Uniform',0,180,5,np);
% for k=1:np
%     model = ManipulatorModel3D(segmentLength,R,r,res, alpha, betas(:,k), ShapeFn);
%     model.fe.plot(model.mesh.nodes);
% end