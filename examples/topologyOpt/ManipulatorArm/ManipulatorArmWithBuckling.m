clear;
close all;

% Computations

% ArmTopOptBucklingFn('Tx', 1.0E7, 0.0E7, 0.0E9, 0.0E7, 0.0E7, 0.0E7); 
% ArmTopOptBucklingFn('Ty', 0.0E7, 1.0E7, 0.0E9, 0.0E7, 0.0E7, 0.0E7);
% ArmTopOptBucklingFn('N',  0.0E7, 0.0E7, 1.0E7, 0.0E7, 0.0E7, 0.0E7);
% ArmTopOptBucklingFn('Mx', 0.0E7, 0.0E7, 0.0E7, 1.0E7, 0.0E7, 0.0E7);
% ArmTopOptBucklingFn('My', 0.0E7, 0.0E7, 0.0E7, 0.0E7, 1.0E7, 0.0E7);
% ArmTopOptBucklingFn('Ms', 0.0E7, 0.0E7, 0.0E7, 0.0E7, 0.0E7, 1.0E7);

% ArmTopOptBucklingFn('My-Ms', 0.0E7, 0.0E7, 0.0E7, 0.0E7, 1.0E7, 1.0E7); 
% ArmTopOptBucklingFn('My-N',  0.0E7, 0.0E7, 1.0E7, 0.0E7, 2.0E7, 0.0E7); 
% ArmTopOptBucklingFn('Ms-N',  0.0E7, 0.0E7, 1.0E7, 0.0E7, 0.0E7, 2.0E7); 

%Plotting results

plotArmTopologies('ArmTopopt_N.mat','Normal force'); savefig('N.fig');
plotArmTopologies('ArmTopopt_Tx.mat','Shear force along x-axis'); savefig('Tx.fig');
plotArmTopologies('ArmTopopt_Ty.mat','Shear force along y-axis'); savefig('Ty.fig');
plotArmTopologies('ArmTopopt_Mx.mat','Bending moment along x-axis'); savefig('Mx.fig');
plotArmTopologies('ArmTopopt_My.mat','Bending moment along y-axis'); savefig('My.fig');
plotArmTopologies('ArmTopopt_Ms.mat','Torsion moment'); savefig('Ms.fig');

plotArmTopologies('ArmTopopt_Ms-N.mat','Mixed mode: torsion and normalforce'); savefig('Ms-N.fig');
plotArmTopologies('ArmTopopt_My-Ms.mat','Mixed mode: y-asis bending and Torsion'); savefig('My-Ms.fig');
plotArmTopologies('ArmTopopt_My-N.mat','Mixed mode: y-asis bending and normal force'); savefig('My-N.fig');
