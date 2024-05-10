clear;
close all;

l=3;         % long edge length    % short edge resolution

b=0.5;
h=2.7;
lc=0.6;
fl=0.5;
hc=0.7;
resb=30;
E=80000;
nu=0.3;

model = CorbelModelLinear(ShapeFunctionL4,b,h,lc,fl,hc,resb,E,nu);

model.plotModel();

randomVariables={RandomVariable("Normal",-2.5,0.1) RandomVariable("Normal",3,0.1)};
transform=IndependentTransformation(randomVariables);
g=performanceFunctionPlus( loadPerformanceFunctionDisp( model ) );


topOpt = StressIntensityTopologyOptimizationVol( 1.2*b/resb, ...
            model.analysis, ... % FEM analysis object
            0.005, ...          % stress intensity treshold for element removal 
            2, ...              % penalty factor
            0.15, ...           % constraint function object
            true ...            % is finite elements uniform
 );

model.setupLoad([-2.5 3]);
topOpt.solve();

% Pdest=[-2.5, 3];

 tuner = ReliabilityTaskTuner(model, topOpt, randomVariables, transform, g, 100000000, 2);
 % tuner.checkTopology(Pdest);
 % [topMin,topMax] = tuner.getBoundsTopologies(0.3,0.6);
 % tuner.fullReliabilityTuning(Pdest,topMin,topMax);



%g.threshold = 0.002;
g.threshold = -0.000812;
% mpps = tuner.checkModality(0.5)
% 
% tuner.tuneMC();
% tuner.plotMCs(["Px" "Py"],'u');

%tuner.fullReliabilityTuning([6, -10]);

tuner.tabPf([-2.5 3],-0.0009,-0.0008,10);

sora2 = SORA('CorbelDispBeta20_2', model,topOpt, randomVariables, g, transform, 2);
sora3 = SORA('CorbelDispBeta20_3', model,topOpt, randomVariables, g, transform, 3);
sora4 = SORA('CorbelDispBeta20_4', model,topOpt, randomVariables, g, transform, 4);
sora5 = SORA('CorbelDispBeta20_5', model,topOpt, randomVariables, g, transform, 5);

% sora2.solveX();
% sora3.solveX();
% sora4.solveX();
% sora5.solveX();


