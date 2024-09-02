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

plotArmTopologies('ArmTopopt_N.mat','Normal force');
plotArmTopologies('ArmTopopt_Tx.mat','Shear force along X direction');
plotArmTopologies('ArmTopopt_Ty.mat','Shear force along Y direction');
plotArmTopologies('ArmTopopt_Mx.mat','Bending moment along X axis');
plotArmTopologies('ArmTopopt_My.mat','Bending moment along Y axis');
plotArmTopologies('ArmTopopt_Ms.mat','Torsion moment');

plotArmTopologies('ArmTopopt_Ms-N.mat','Mixed mode Torsion and normalforce');
plotArmTopologies('ArmTopopt_My-Ms.mat','Mixed mode Y asis bending and Torsion');
plotArmTopologies('ArmTopopt_My-N.mat','Mixed mode Y asis bending and normal force');
