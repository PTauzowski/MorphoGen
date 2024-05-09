clear;
close all;


resb=20;
E=80000;
nu=0.3;

model = CorbelModelLinear3D(ShapeFunctionL8,resb,E,nu);

model.plotModel();
view(45,45);

randomVariables={RandomVariable("Normal",-2.5,0.25)};
transform=IndependentTransformation(randomVariables);
g=performanceFunctionPlus( loadPerformanceFunctionDisp( model ) );


topOpt = StressIntensityTopologyOptimizationVol( 1.2*0.05/resb, ...
            model.analysis, ... % FEM analysis object
            0.005, ...          % stress intensity treshold for element removal 
            2, ...              % penalty factor
            0.4, ...           % constraint function object
            true ...            % is finite elements uniform
 );

model.setupLoad([-2.5 3]);
topOpt.solve();

Pdest=[-2.5, 3];

 tuner = ReliabilityTaskTuner(model, topOpt, randomVariables, transform, g, 100000000, 2);
 % tuner.checkTopology(Pdest);
 % [topMin,topMax] = tuner.getBoundsTopologies(0.3,0.6);
 % tuner.fullReliabilityTuning(Pdest,topMin,topMax);


g.threshold = 0.000;
%g.threshold = 0.00;


tuner.tuneMC();
tuner.plotMCs(["Px" "Py"],'u');

%tuner.fullReliabilityTuning([6, -10]);

sora2 = SORA('LShapeDispBeta30_2', model,topOpt, randomVariables, g, transform, 2);
sora3 = SORA('LShapeDispBeta30_3', model,topOpt, randomVariables, g, transform, 3);
sora4 = SORA('LShapeDispBeta30_4', model,topOpt, randomVariables, g, transform, 4);
sora5 = SORA('LShapeDispBeta30_5', model,topOpt, randomVariables, g, transform, 5);

 % sora2.solveX();
 % sora3.solveX();
 % sora4.solveX();
 % sora5.solveX();


