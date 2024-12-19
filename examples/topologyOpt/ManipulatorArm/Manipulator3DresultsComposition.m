clear;
close all;

load("ManipulatorBucklingAndTorsion.mat");

% plotArmConfigurationHMextended("maximal Huber-Mises stress configuration 1, \sigma_{HM}=" + num2str(vSort(nSamples-0)),'MaxHM_configuration1.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSort(nSamples-0),:), ShapeFn, 0.0);
% plotArmConfigurationHMextended("maximal Huber-Mises stress configuration 2, \sigma_{HM}=" + num2str(vSort(nSamples-1)),'MaxHM_configuration2.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSort(nSamples-1),:), ShapeFn, 0.0);
% plotArmConfigurationHMextended("maximal Huber-Mises stress configuration 3, \sigma_{HM}=" + num2str(vSort(nSamples-2)),'MaxHM_configuration3.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSort(nSamples-2),:), ShapeFn, 0.0);
% 
% plotArmConfigurationHMextended("minimal Huber-Mises stress configuration 1, \sigma_{HM}=" + num2str(vSort(1)),'MinHM_configuration1.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSort(1),:), ShapeFn, 0.0);
% plotArmConfigurationHMextended("minimal Huber-Mises stress configuration 2, \sigma_{HM}=" + num2str(vSort(2)),'MinHM_configuration2.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSort(2),:), ShapeFn, 0.0);
% plotArmConfigurationHMextended("minimal Huber-Mises stress configuration 3, \sigma_{HM}=" + num2str(vSort(3)),'MinHM_configuration3.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSort(3),:), ShapeFn, 0.0);
% 
% plotArmConfigurationHMextended("maximal torsion moment configuration 1, M_{s}=" + num2str(vSortMs1(nSamples-0)),'MaxMs_configuration1.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSortMs1(nSamples-0),:), ShapeFn, 0.0);
% plotArmConfigurationHMextended("maximal torsion moment configuration 2, M_{s}=" + num2str(vSortMs1(nSamples-1)),'MaxMs_configuration2.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSortMs1(nSamples-1),:), ShapeFn, 0.0);
% plotArmConfigurationHMextended("maximal torsion moment configuration 3, M_{s}=" + num2str(vSortMs1(nSamples-2)),'MaxMs_configuration3.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSortMs1(nSamples-2),:), ShapeFn, 0.0);
% 
% plotArmConfigurationHMextended("minimal torsion moment configuration 1, M_{s}=" + num2str(vSortMs1(1)),'MinMs_configuration1.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSortMs1(1),:), ShapeFn, 0.0);
% plotArmConfigurationHMextended("minimal torsion moment configuration 2, M_{s}=" + num2str(vSortMs1(2)),'MinMs_configuration2.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSortMs1(2),:), ShapeFn, 0.0);
% plotArmConfigurationHMextended("minimal torsion moment configuration 3, M_{s}=" + num2str(vSortMs1(3)),'MinMs_configuration3.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSortMs1(3),:), ShapeFn, 0.0);


%plotArmConfigurationHMextended("maximal Huber-Mises stress configuration 1, \sigma_{HM}=" + num2str(vSort(nSamples-0)),'BendingConfiguration.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, 1, alpha, samples(iSort(nSamples-0),:), ShapeFn, 0.0);
%plotArmConfigurationHMextended("maximal torsion moment configuration 1, M_{s}=" + num2str(vSortMs1(nSamples-0)),'TorsionConfiguration.pdf',E,nu,segmentLength,R,r,res, modelMax.hal, fSegmentNelems, 2, alpha, samples(iSortMs1(nSamples-0),:), ShapeFn, 0.0);

plotArmTopOptConfigProjections("BendingTopology","Bending topology",Rfilter, modelMaxMs.analysis, modelMax.halfSegmentNelems, 1, xopt_bending_buckling, cutTreshold, penal, false);
%plotArmTopOptConfigProjections("TorsionTopology","Torsion topology",Rfilter, modelMaxMs.analysis, modelMax.halfSegmentNelems, 2, xopt_torsion_buckling, cutTreshold, penal, false);

