clear;
close all;

l=3;         % long edge length    % short edge resolution

b=0.5;
h=2.7;
lc=0.6;
fl=0.5;
hc=0.7;
resb=20;
E1=80000;
E2=80000;
E3=80000;
P1=-2.5;
P2=3;
nu=0.3;

model = CorbelModelMultiMat(ShapeFunctionL4,b,h,lc,fl,hc,resb,E1,E2,E3,nu);
model.plotModel();


randomVariables={RandomVariable("Normal",E1,0.1*E1) RandomVariable("Normal",E2,0.1*E2) RandomVariable("Normal",E3,0.1*E3),RandomVariable("Normal",P1,0.01*abs(P1)) RandomVariable("Normal",P2,0.01*abs(P2))};
transform=IndependentTransformation(randomVariables);
g=performanceFunctionPlus( loadPerformanceFunctionDisp( model, length(randomVariables) ) );
%g=performanceFunctionPlus( penalizedStressPerformanceFunction( model,length(randomVariables),6) );




topOpt = StressIntensityTopologyOptimizationVol( 1.2*b/resb, ...
            model.analysis, ... % FEM analysis object
            0.005, ...          % stress intensity treshold for element removal 
            2, ...              % penalty factor
            0.2, ...           % constraint function object
            true ...            % is finite elements uniform
 );

model.setupLoad([ randomVariables{4}.mean randomVariables{5}.mean ] );
% model.solveWeighted();
% 
% model.plotModel();
% model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
% for k=1:size(model.fe,2)
%     model.fe{k}.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);
% end

 % 
 % topOpt.solve();
 % title("Starting topology");
 % figure;

% Pdest=[-2.5, 3];

 tuner = ReliabilityTaskTuner(model, topOpt, randomVariables, transform, g, 2000, 2);
 %tuner = ReliabilityTaskTuner(model, topOpt, randomVariables, transform, g, 2000, 2);
 % tuner.checkTopology(Pdest);
 % [topMin,topMax] = tuner.getBoundsTopologies(0.3,0.6);
 % tuner.fullReliabilityTuning(Pdest,topMin,topMax);


% g.threshold = 0.028;
g.threshold = 30;
%g.threshold = 6.63;
% mpps = tuner.checkModality(0.5)


 % tuner.tuneMC();
 % tuner.plotMCs(["E1" "E2" "E3" "Px" "Py"],'u');
 tuner.tabPf()

% sora2 = SORAold('CorbelStress20_2', model,topOpt, randomVariables, g, transform, 2);
% sora3 = SORAold('CorbelStress20_3', model,topOpt, randomVariables, g, transform, 3);
% sora4 = SORAold('CorbelStress20_4', model,topOpt, randomVariables, g, transform, 4);
% sora5 = SORAold('CorbelStress20_5', model,topOpt, randomVariables, g, transform, 5);

sora2 = SORAold('CorbelDisp50_2', model,topOpt, randomVariables, g, transform, 2);
sora3 = SORAold('CorbelDisp50_3', model,topOpt, randomVariables, g, transform, 3);
sora4 = SORAold('CorbelDisp50_4', model,topOpt, randomVariables, g, transform, 4);
sora5 = SORAold('CorbelDisp50_5', model,topOpt, randomVariables, g, transform, 5);


% sora5.solveX(0.38);
% sora2.solveX(0.25);
% sora3.solveX(0.34);
% sora4.solveX(0.36);


