clear;
close all;
height=12;
E=210000;
nu=0.3;
alphaT=1 ;
dT=-20;

 % x(1) - gan, alGan proportion,  x(2) - relNotchDepth ,    x(3) - notchWidth
 % 0.2358    0.4383    5.0757 - test
lbg=[0.2 0.1 1];
ubg=[0.8 0.9 8];
x0g=(lbg+ubg)/2;

g = chocolatePerformanceFunction(height,210000,0.3,alphaT,dT,x0g);
ntv=g.model.nTempVars;
lbt=zeros(1,ntv);
ubt=ones(1,ntv);
x0t=ones(1,ntv)*0.5;

lb=[lbg lbt];
ub=[ubg ubt];
x0=[x0g x0t];

fn_g = @(x)( g.computeValue(x) );
%g.fullFactorialBoundsPlot(lb,ub);
N=8000;
randomVariables=cell(1,size(x0,2));
for k=1:size(x0,2)
    randomVariables{k}=RandomVariable("Uniform",lb(k),ub(k));
end
transform=IndependentTransformation(randomVariables);
mc= MonteCarlo(randomVariables,g,N);
%x = mc.generateRandomSapmles(N);
tic
%res_mc = mc.solve();
x0(1:3)=[0.3333    0.8    4  ];
x0(4:6)=[1    1    1  ];
xopt = fmincon(fn_g,x0,[],[],[],[],lb,ub)
%xopt = fmincon(fn_g,x0)
toc

save("chocolateOptiThermalBestX0.mat");

% [v, i]=min(mc.r)
% g.evaluateValue(mc.x(i,:))
%x0(1:3)=[0.7830    0.8724    2.7589  ];

%g.evaluateValue(x0);
g.evaluateValue(xopt);
%g.createModel([0.2358    0.4383    5.0757]);
%g.evaluateValue2();



