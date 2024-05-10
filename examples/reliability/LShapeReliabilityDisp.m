clear;
close all;

l=3;         % long edge length
res = 40;    % short edge resolution

model = LShapeModelLinear(  ShapeFunctionL4,...
                             l, ...         % long edge length
                            res, ...         % short edge resolution
                            210000, ...     % Young Modulus
                            0.3, ...        % Poinsson's ratio
                            [l 0.4*l] ...   % Force location,xp
     );

model.setResultNode([l l*0.2]);
model.plotModel();

% [5.8 -7.7
%randomVariables={RandomVariable("Lognormal",-1,1) RandomVariable("Normal",-6,2)};
%randomVariables={RandomVariable("Normal",-1,4) RandomVariable("Normal",-4,1)};
randomVariables={RandomVariable("Normal",7,1) RandomVariable("Normal",-4,1)};
%randomVariables={RandomVariable("Normal",0,2) RandomVariable("Normal",-10,1)};
transform=IndependentTransformation(randomVariables);
g=performanceFunctionPlus( loadPerformanceFunctionDisp( model ) );


topOpt = StressIntensityTopologyOptimizationVol( 1.2*l/res, ...
            model.analysis, ... % FEM analysis object
            0.005, ...          % stress intensity treshold for element removal 
            2, ...              % penalty factor
            0.2, ...           % constraint function object
            true ...            % is finite elements uniform
 );

model.setupLoad([7 -4]);
topOpt.solve();

Pdest=[6, -10];

 tuner = ReliabilityTaskTuner(model, topOpt, randomVariables, transform, g, 100000000, 2);
 % tuner.checkTopology(Pdest);
 % [topMin,topMax] = tuner.getBoundsTopologies(0.3,0.6);
 % tuner.fullReliabilityTuning(Pdest,topMin,topMax);



g.threshold = 0.0055;
%g.threshold = 0.00;


 %mpps = tuner.checkModality(0.5)

%tuner.tuneMC();
% tuner.plotMCs(["Px" "Py"],'u');

%tuner.fullReliabilityTuning([6, -10]);

sora2 = SORA('LShapeDispBeta40c_2', model,topOpt, randomVariables, g, transform, 2);
sora3 = SORA('LShapeDispBeta40c_3', model,topOpt, randomVariables, g, transform, 3);
sora4 = SORA('LShapeDispBeta40c_4', model,topOpt, randomVariables, g, transform, 4);
sora5 = SORA('LShapeDispBeta40c_5', model,topOpt, randomVariables, g, transform, 5);

 sora2.solveX();
 sora3.solveX();
 sora4.solveX();
 sora5.solveX();


