clear;
close all;

l=3; % long edge length

model = LShapeModelLinear(  ShapeFunctionL4,...
                             l, ...         % long edge length
                            40, ...         % short edge resolution
                            210000, ...     % Young Modulus
                            0.3, ...        % Poinsson's ratio
                            [l 0.4*l] ...   % Force location,xp
     );

model.setResultNode([l l*0.4]);
model.plotModel();

randomVariables={RandomVariable("Normal",0,3) RandomVariable("Normal",10,1)};
transform=IndependentTransformation(randomVariables);
g=loadPerformanceFunctionDisp(model);


topOpt = StressIntensityTopologyOptimization( Rfilter, ...
            model.analysis, ... % FEM analysis object
            0.005, ...          % stress intensity treshold for element removal 
            2, ...              % penalty factor
            g, ...              % constraint function object
            true ...            % is finite elements uniform
 );

[objF, xopt]  = topOpt.solve();

sora2 = SORA('LShapeDispBeta20_2', model,topOpt, randomVariables, g, transform, 2);
sora3 = SORA('LShapeDispBeta20_3', model,topOpt, randomVariables, g, transform, 3);
sora4 = SORA('LShapeDispBeta20_4', model,topOpt, randomVariables, g, transform, 4);
sora5 = SORA('LShapeDispBeta20_5', model,topOpt, randomVariables, g, transform, 5);

sora2.solveX();
sora3.solveX();
sora4.solveX();
sora5.solveX();


