clear;
close all;

l=3;         % long edge length
res = 40;    % short edge resolution

model = LShapeModelLinear(  ShapeFunctionL4,...
                             l, ...         % long edge length
                            res, ...         % short edge resolution
                            210000, ...     % Young Modulus
                            0.3, ...        % Poinsson's ratio
                            [l 0.2*l] ...   % Force location,xp
     );

model.setResultNode([l l*0.2]);
model.plotModel();

randomVariables={RandomVariable("Normal",5,2) RandomVariable("Normal",-10,2)};
transform=IndependentTransformation(randomVariables);
g=loadPerformanceFunctionDisp(model);


topOpt = StressIntensityTopologyOptimizationVol( 1.2*l/res, ...
            model.analysis, ... % FEM analysis object
            0.005, ...          % stress intensity treshold for element removal 
            2, ...              % penalty factor
            0.4, ...              % constraint function object
            true ...            % is finite elements uniform
 );


% tuner = ReliabilityTaskTuner(model, topOpt, randomVariables, transform, g, 1000000, 2);
% tuner.tuneMC();
% tuner.plotMCs(["Px" "Py"],'u');

sora2 = SORA('LShapeDispBeta20_2', model,topOpt, randomVariables, g, transform, 2);
sora3 = SORA('LShapeDispBeta20_3', model,topOpt, randomVariables, g, transform, 3);
sora4 = SORA('LShapeDispBeta20_4', model,topOpt, randomVariables, g, transform, 4);
sora5 = SORA('LShapeDispBeta20_5', model,topOpt, randomVariables, g, transform, 5);

sora2.solveX();
sora3.solveX();
sora4.solveX();
sora5.solveX();


