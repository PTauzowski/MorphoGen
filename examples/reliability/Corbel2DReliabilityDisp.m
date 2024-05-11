clear;
close all;

l=3;         % long edge length    % short edge resolution

b=0.5;
h=2.7;
lc=0.6;
fl=0.5;
hc=0.7;
resb=50;
E=80000;
nu=0.3;

model = CorbelModelLinear(ShapeFunctionL4,b,h,lc,fl,hc,resb,E,nu);

model.plotModel();

randomVariables={RandomVariable("Normal",-2.5,0.25) RandomVariable("Normal",3,0.3)};
transform=IndependentTransformation(randomVariables);
% g=performanceFunctionPlus( loadPerformanceFunctionDisp( model ) );
g=performanceFunctionMinus( penalizedStressPerformanceFunction(model,6) );


topOpt = StressIntensityTopologyOptimizationVol( 1.2*b/resb, ...
            model.analysis, ... % FEM analysis object
            0.005, ...          % stress intensity treshold for element removal 
            2, ...              % penalty factor
            0.15, ...           % constraint function object
            true ...            % is finite elements uniform
 );
 % 
 % model.setupLoad([-2.5 3]);
 % topOpt.solve();

% Pdest=[-2.5, 3];

 %tuner = ReliabilityTaskTuner(model, topOpt, randomVariables, transform, g, 100000000, 2);
 tuner = ReliabilityTaskTuner(model, topOpt, randomVariables, transform, g, 2000, 2);
 % tuner.checkTopology(Pdest);
 % [topMin,topMax] = tuner.getBoundsTopologies(0.3,0.6);
 % tuner.fullReliabilityTuning(Pdest,topMin,topMax);


%g.threshold = 0.002;
g.threshold = 7.5;
% mpps = tuner.checkModality(0.5)


% tuner.tuneMC();
% tuner.plotMCs(["Px" "Py"],'u');


sora2 = SORA('CorbelStress50_2', model,topOpt, randomVariables, g, transform, 2);
sora3 = SORA('CorbelStress50_3', model,topOpt, randomVariables, g, transform, 3);
sora4 = SORA('CorbelStress50_4', model,topOpt, randomVariables, g, transform, 4);
sora5 = SORA('CorbelStress50_5', model,topOpt, randomVariables, g, transform, 5);

% tuner.tabPf([-2.5 3],0,0.06,10);

sora2.solveX();
sora5.solveX();
sora3.solveX();
sora4.solveX();



