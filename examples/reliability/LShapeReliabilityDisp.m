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

% [5.8 -7.7
%randomVariables={RandomVariable("Lognormal",-1,1) RandomVariable("Normal",-6,2)};
%randomVariables={RandomVariable("Normal",-1,4) RandomVariable("Normal",-4,1)};
randomVariables={RandomVariable("Normal",3,3) RandomVariable("Normal",-10,1)};
transform=IndependentTransformation(randomVariables);
g=loadPerformanceFunctionDisp(model);


topOpt = StressIntensityTopologyOptimizationVol( 1.2*l/res, ...
            model.analysis, ... % FEM analysis object
            0.005, ...          % stress intensity treshold for element removal 
            2, ...              % penalty factor
            0.2, ...           % constraint function object
            true ...            % is finite elements uniform
 );

%model.setupLoad([6, -10]);
%topOpt.solve();


% tuner = ReliabilityTaskTuner(model, topOpt, randomVariables, transform, g, 1000000, 2);
% tuner.tuneMC();
% tuner.plotMCs(["Px" "Py"],'u');

sora2 = SORA('LShapeDispBeta20_2', model,topOpt, randomVariables, g, transform, 2);
sora3 = SORA('LShapeDispBeta20_3', model,topOpt, randomVariables, g, transform, 3);
sora4 = SORA('LShapeDispBeta20_4', model,topOpt, randomVariables, g, transform, 4);
sora5 = SORA('LShapeDispBeta20_5', model,topOpt, randomVariables, g, transform, 5);

% %sora0.solveX();
% sora2.solveX();
% sora3.solveX();
% sora4.solveX();
sora5.solveX();


