clear;
close all;

% Cantilever topology optimization elastic task.

% Resolution of shortest (vertical) edge
res = 40;

% height of the cantilever
h = 1;

% Aspect ratio length/height
aspect=6;

% length of the cantilever
l=aspect*h;

% Filtering radius
Rfilter = 2*h/res;

%Removal intensity threshold
cutTreshold = 0.005;

%penalty factor
penal = 3;

% Type of shape function to be used (here: four node Langrange)
sfL4 = ShapeFunctionL4;

% Creating FE mesh object
mesh = Mesh();

% Generating rectangular mesh ( aspect*h x h )
mesh.addRectMesh2D(0, 0, l, h, aspect*res, res, sfL4.pattern);

% Create plane stress finite element object
fe=PlaneStressElem( sfL4, mesh.elems );

% Create isotropic material object
material = PlaneStressMaterial('mat1');
material.setElasticIzo(205E9, 0.3);
material.setMassIzoMatrix(7850);
% Assigning material to finite element
fe.setMaterial( material );

% Ke = fe.computeStifnessMatrix(mesh.nodes,1,1);
% Me = fe.computeMassMatrix(mesh.nodes,1,1);
% 
%  [q, l] = eigs( Ke(:,:,1),Me(:,:,1), 8, 'smallestabs');

% Creating linear elastic finite element analysis object with weighted matrix feature, weighted by element density.
analysisLinear      = LinearElasticityWeighted( fe, mesh, false );
analysisSecondOrder = SecondOrderDynamicElasticityWeighted(fe, mesh, 0, false);

% Creating node selector object to select fixed edge (left)
loadEdgeSelector = Selector( @(x)( abs(x(:,2) -h ) < 0.0005 ) );
fixedEdgeSelector = Selector( @(x)( abs(x(:,1)) < 0.0005 ) );
fixedEdgeSelector2 = Selector( @(x)( abs(x(:,1)-l) < 0.0005 ) );

% Fixing structure according to above defined node selector object
analysisLinear.fixNodes( fixedEdgeSelector, ["ux" "uy"] );
analysisSecondOrder.fixNodes( fixedEdgeSelector2, ["ux"] );
analysisSecondOrder.fixNodes( fixedEdgeSelector, ["ux" "uy"] );

% Creating load vector with one node loaded at the middle of right edge
P=-2.0E9; %100;
%P=-1.5E8; %150;

hp=h/2;
analysisLinear.loadClosestNode([aspect*h, hp ], ["ux" "uy"], [0 P] );
%analysisSecondOrder.loadClosestNode([aspect*h, hp ], ["ux" "uy"], [0 P] );
analysisSecondOrder.elementLoadLineIntegral( "global", loadEdgeSelector,  ["ux" "uy"], @(x)( x*0 + [0 P/l] ));


const_rows=3;
ncel=round(const_rows*res*l);
const_elems=[1:ncel size(mesh.elems,1):-1:size(mesh.elems,1)-ncel ];

nEigenForms=10;
stability = LinearStability( analysisLinear.felems, mesh);
stability.Pnodal = analysisLinear.Pnodal;
stability.Pfem = analysisLinear.Pfem;
stability.supports = analysisLinear.supports;
% stability.solve( nEigenForms);
% lambdas = diag(stability.lambdas);
% for k=1:min(10,nEigenForms)
%     subplot(5, 2, k);
%     stability.setForm(k);
%     fe.plotWired(mesh.nodes,stability.qnodal,0.2);
%     axis on, xlabel('x-axis'), ylabel('y-axis'), view(3)
%     lambda_str = sprintf('%.4g', lambdas(k));
%     title(['Form:' num2str(k), ' \lambda=' lambda_str]);
% end

nEigenForms=10;
vibrations = LinearNaturalVibration( analysisLinear.felems, mesh );
vibrations.Pnodal = analysisLinear.Pnodal;
vibrations.Pfem = analysisLinear.Pfem;
vibrations.supports = analysisLinear.supports;
% vibrations.solve( nEigenForms);
% omegas = diag(vibrations.omegas)
% for k=1:min(10,nEigenForms)
%     %subplot(5, 2, k);
%     figure
%     vibrations.setForm(k);
%     fe.plotWithSettings(mesh.nodes,"deformed",vibrations.qnodal,0.2);
%     %axis on, xlabel('x-axis'), ylabel('y-axis'), view(3)
%     omega_str = sprintf('%.4g', omegas(k));
%     title(['Form:' num2str(k), ' \omega=' omega_str]);
% end

analysisWithBuckling = SecondOrderDynamicElasticityWeighted( fe, mesh, 0.30, false );
analysisWithBuckling.Pnodal=stability.Pnodal;
analysisWithBuckling.Pfem=stability.Pfem;
analysisWithBuckling.supports=stability.supports;

harmonicVivrations = ElasticHarmonicVibrations(fe, mesh, 2*pi*20, true);
harmonicVivrations.Pnodal=analysisSecondOrder.Pnodal;
harmonicVivrations.Pfem=analysisSecondOrder.Pfem;
harmonicVivrations.supports=analysisSecondOrder.supports;


% figure;
% tic
% topOptLinear = StressIntensityTopologyOptimizationVol( Rfilter, analysisLinear, cutTreshold, penal, 0.4, true );
% [objF, xopt]  = topOptLinear.solve();
% toc

figure;
tic
topOptHarmonicVivrations = StressIntensityTopologyOptimizationVol( Rfilter, harmonicVivrations, cutTreshold, penal, 0.4, true );
topOptHarmonicVivrations.setConstElems(const_elems);
[objF, xopt]  = topOptHarmonicVivrations.solve();
toc

% figure;
% tic
% topOptSecondOrder = StressIntensityTopologyOptimizationDynamicBuckling( Rfilter, analysisSecondOrder, cutTreshold, penal, 0.30, true );
% topOptSecondOrder.setConstElems(const_elems);
% [objF, xopt]  = topOptSecondOrder.solve();
% toc

% figure;
% tic
% topOptBuckling = StressIntensityTopologyOptimizationBuckling( Rfilter, analysisWithBuckling, cutTreshold, penal, 0.38, true );
% [objF, xopt]  = topOptBuckling.solve();
% toc

save('RampZVibrations.mat');

figure, hold on
p1=plot(topOptSecondOrder.plVol,topOptSecondOrder.plLambda,'b','LineWidth', 3);
set(gca, 'XDir', 'reverse');
title('Critical force coefficient evolution without buckling');
xlabel('Volume fracion [%]');
ylabel('Critical force coefficient [%]');
%xlim([37 57]);
set(gca, 'FontSize', 18)

figure, hold on
p2=plot(topOptSecondOrder.plVol,topOptSecondOrder.plLambda,'r','LineWidth', 3);
set(gca, 'XDir', 'reverse');
title('Critical force coefficient evolution with buckling');
xlabel('Volume fracion [%]');
ylabel('Critical force coefficient [%]');
%xlim([37 57]);
set(gca, 'FontSize', 24)


legend_modes={'Mode 1', 'Mode 2', 'Mode 3', 'Mode 4', 'Mode 5'};
for l=1:5
    figure, hold on
    p2=plot(topOptSecondOrder.plVol,topOptSecondOrder.plOmegas(l,:)','LineWidth', 3);
    set(gca, 'XDir', 'reverse');
    title(['Mode '   num2str(l)  ' evolution']);
    %% 
    xlabel('Volume fracion [%]');
    ylabel('Frequency [Hz]');
    %xlim([37 57]);
    set(gca, 'FontSize', 18)
    saveas(gcf,['Mode_'  num2str(l)  '.png'])
    savefig(gcf,['Mode_'  num2str(l)  '.fig'])
end


fontsize=16;
iters=[1 100 250 350  450 size(topOptSecondOrder.allx,2)-2];
for i=iters(1):size(iters,2)
    figure;
    subplot(2, 1, 1);
    fe.plotWithSettings(mesh.nodes,"deformed",analysisSecondOrder.fromFEMVector( analysisSecondOrder.modes1(:,i) ),0.2,"elem nums",topOptSecondOrder.allx(:,i)>=0.5);
    title(['Mode 1, vol_{fr}=' num2str(topOptSecondOrder.plVol(i)) ', Frq. =' num2str(topOptSecondOrder.plOmegas(1,i),4) ' [Hz]'  ', iter:' num2str(i)]);
    set(gca, 'FontSize', fontsize)
    
    subplot(2, 1, 2);
    fe.plotWithSettings(mesh.nodes,"deformed",analysisSecondOrder.fromFEMVector( analysisSecondOrder.modes2(:,i) ),0.2,"elem nums",topOptSecondOrder.allx(:,i)>=0.5);
    title(['Mode 2, vol_{fr}=' num2str(topOptSecondOrder.plVol(i)) ', Frq. =' num2str(topOptSecondOrder.plOmegas(2,i),4) ' [Hz]' ', iter' num2str(i)]);
    set(gca, 'FontSize', fontsize)
    saveas(gcf,['frame_' num2str(i) '.pdf'])
    savefig(gcf,['frame_' num2str(i) '.fig'])
end

% figure;
% U = normalize(analysisSecondOrder.modes2, 1); 
% C = abs(U' * U);
% imagesc(C);
% colormap(flipud(hot));
% colorbar;
% axis equal tight;
% xlabel('Mode number');
% ylabel('Mode number');
% title('Correlation matrix of eigenmode 1');
% 
% figure;
% U = normalize(analysisSecondOrder.modes2, 1); 
% C = abs(U' * U);
% imagesc(C);
% colormap(flipud(hot));
% colorbar;
% axis equal tight;
% xlabel('Mode number');
% ylabel('Mode number');
% title('Correlation matrix of eigenmode 2');
% 
% figure;
% U = normalize(analysisSecondOrder.modes3, 1); 
% C = abs(U' * U);
% imagesc(C);
% colormap(flipud(hot));
% colorbar;
% axis equal tight;
% xlabel('Mode number');
% ylabel('Mode number');
% title('Correlation matrix of eigenmode 3');

pl_mode1_cor = analysisSecondOrder.getModesCorrelation(analysisSecondOrder.modes1);
pl_mode2_cor = analysisSecondOrder.getModesCorrelation(analysisSecondOrder.modes2);
pl_mode3_cor = analysisSecondOrder.getModesCorrelation(analysisSecondOrder.modes3);
pl_mode4_cor = analysisSecondOrder.getModesCorrelation(analysisSecondOrder.modes4);
pl_mode5_cor = analysisSecondOrder.getModesCorrelation(analysisSecondOrder.modes5);

pl_modes_cor = [analysisSecondOrder.getModesCorrelation(analysisSecondOrder.modes1); analysisSecondOrder.getModesCorrelation(analysisSecondOrder.modes2); analysisSecondOrder.getModesCorrelation(analysisSecondOrder.modes3); analysisSecondOrder.getModesCorrelation(analysisSecondOrder.modes4); analysisSecondOrder.getModesCorrelation(analysisSecondOrder.modes5)];


for k=1:5
    kstr=num2str(k);
    figure, hold on
    p2=plot(topOptSecondOrder.plVol(1:end-1), pl_modes_cor(k,1:end)','LineWidth', 2);
    set(gca, 'XDir', 'reverse');
    title(['Correlations of mode ' kstr]);
    xlabel('Volume fracion [%]');
    ylabel('Correlation');
    %xlim([37 57]);
    set(gca, 'FontSize', 18)
    saveas(gcf,['correlation_' kstr '.png'])
    savefig(gcf,['correlation_' kstr '.fig'])
end

for k=1:size(pl_mode1_cor,2)
    if pl_mode1_cor(k)<0.3
        figure;

        subplot(2, 1, 1);
        fe.plotWithSettings(mesh.nodes,"deformed",analysisSecondOrder.fromFEMVector( analysisSecondOrder.modes1(:,k) ),0.1,"elem nums",topOptSecondOrder.allx(:,k)>=0.5,"nodes",false);
        title(['Mode 1, vol_{fr}=' num2str(topOptSecondOrder.plVol(k)) ', Frq. =' num2str(topOptSecondOrder.plOmegas(1,k-1),4) ' [Hz]'], [ 'MAC=' num2str(pl_mode1_cor(k),3) ', iter:' num2str(k)]);
        set(gca, 'FontSize', fontsize)
        
        subplot(2, 1, 2);
        fe.plotWithSettings(mesh.nodes,"deformed",analysisSecondOrder.fromFEMVector( analysisSecondOrder.modes1(:,k+1) ),0.1,"elem nums",topOptSecondOrder.allx(:,k+1)>=0.5,"nodes",false);
        title(['Mode 1, vol_{fr}=' num2str(topOptSecondOrder.plVol(k+1)) ', Frq. =' num2str(topOptSecondOrder.plOmegas(1,k+1),4) ' [Hz]' ', iter' num2str(k+1)]);
        set(gca, 'FontSize', fontsize)

        saveas(gcf,['frame_correlation_mode_1_' num2str(k) '.pdf'])
        savefig(gcf,['frame_correlation_mode_1_' num2str(k) '.fig'])
    end
end

for k=1:size(pl_mode2_cor,2)
    if pl_mode2_cor(k)<0.3
        figure;
        subplot(2, 1, 1);
        fe.plotWithSettings(mesh.nodes,"deformed",analysisSecondOrder.fromFEMVector( analysisSecondOrder.modes2(:,k) ),0.1,"elem nums",topOptSecondOrder.allx(:,k)>=0.5,"nodes",false);
        title(['Mode 2, vol_{fr}=' num2str(topOptSecondOrder.plVol(k)) ', Frq. =' num2str(topOptSecondOrder.plOmegas(2,k),4) ' [Hz]' ], [ 'MAC=' num2str(pl_mode2_cor(k),3) ', iter:' num2str(k)]);
        set(gca, 'FontSize', fontsize)
        
        subplot(2, 1, 2);
        fe.plotWithSettings(mesh.nodes,"deformed",analysisSecondOrder.fromFEMVector( analysisSecondOrder.modes2(:,k+1) ),0.1,"elem nums",topOptSecondOrder.allx(:,k+1)>=0.5,"nodes",false);
        title(['Mode 2, vol_{fr}=' num2str(topOptSecondOrder.plVol(k+1)) ', Frq. =' num2str(topOptSecondOrder.plOmegas(2,k+1),4) ' [Hz]' ', iter' num2str(k+1)]);
        set(gca, 'FontSize', fontsize)
        saveas(gcf,['frame_correlation_mode_2_' num2str(k) '.pdf'])
        savefig(gcf,['frame_correlation_mode_2_' num2str(k) '.fig'])
    end
end

for k=1:size(pl_mode3_cor,2)
    if pl_mode3_cor(k)<0.3
        figure;
        subplot(2, 1, 1);
        fe.plotWithSettings(mesh.nodes,"deformed",analysisSecondOrder.fromFEMVector( analysisSecondOrder.modes3(:,k) ),0.1,"elem nums",topOptSecondOrder.allx(:,k)>=0.5,"nodes",false);
        title(['Mode 3, vol_{fr}=' num2str(topOptSecondOrder.plVol(k)) ', Frq. =' num2str(topOptSecondOrder.plOmegas(3,k),4) ' [Hz]' ],[ 'MAC=' num2str(pl_mode3_cor(k),3) ', iter:' num2str(k)]);
        set(gca, 'FontSize', fontsize)
        
        subplot(2, 1, 2);
        fe.plotWithSettings(mesh.nodes,"deformed",analysisSecondOrder.fromFEMVector( analysisSecondOrder.modes3(:,k+1) ),0.1,"elem nums",topOptSecondOrder.allx(:,k+1)>=0.5,"nodes",false);
        title(['Mode 3, vol_{fr}=' num2str(topOptSecondOrder.plVol(k+1)) ', Frq. =' num2str(topOptSecondOrder.plOmegas(3,k+1),4) ' [Hz]' ', iter' num2str(k+1)]);
        set(gca, 'FontSize', fontsize)
        saveas(gcf,['frame_correlation_mode_3_' num2str(k) '.pdf'])
        savefig(gcf,['frame_correlation_mode_3_' num2str(k) '.fig'])
    end
end

for k=1:size(pl_mode4_cor,2)
    if pl_mode4_cor(k)<0.3
        figure;
        subplot(2, 1, 1);
        fe.plotWithSettings(mesh.nodes,"deformed",analysisSecondOrder.fromFEMVector( analysisSecondOrder.modes4(:,k) ),0.1,"elem nums",topOptSecondOrder.allx(:,k)>=0.5,"nodes",false);
        title(['Mode 4, vol_{fr}=' num2str(topOptSecondOrder.plVol(k)) ', Frq. =' num2str(topOptSecondOrder.plOmegas(4,k),4) ' [Hz]'],[ 'MAC=' num2str(pl_mode4_cor(k),3) ', iter:' num2str(k)]);
        set(gca, 'FontSize', fontsize)
        
        subplot(2, 1, 2);
        fe.plotWithSettings(mesh.nodes,"deformed",analysisSecondOrder.fromFEMVector( analysisSecondOrder.modes3(:,k+1) ),0.1,"elem nums",topOptSecondOrder.allx(:,k+1)>=0.5,"nodes",false);
        title(['Mode 4, vol_{fr}=' num2str(topOptSecondOrder.plVol(k+1)) ', Frq. =' num2str(topOptSecondOrder.plOmegas(4,k+1),4) ' [Hz]' ', iter' num2str(k+1)]);
        set(gca, 'FontSize', fontsize)
        saveas(gcf,['frame_correlation_mode_4_' num2str(k) '.pdf'])
        savefig(gcf,['frame_correlation_mode_4_' num2str(k) '.fig'])
    end
end

for k=1:size(pl_mode5_cor,2)
    if pl_mode5_cor(k)<0.3
        figure;
        subplot(2, 1, 1);
        fe.plotWithSettings(mesh.nodes,"deformed",analysisSecondOrder.fromFEMVector( analysisSecondOrder.modes5(:,k) ),0.1,"elem nums",topOptSecondOrder.allx(:,k)>=0.5,"nodes",false);
        title(['Mode 5, vol_{fr}=' num2str(topOptSecondOrder.plVol(k)) ', Frq. =' num2str(topOptSecondOrder.plOmegas(5,k),4) ' [Hz]' ],[ 'MAC=' num2str(pl_mode3_cor(k),3) ', iter:' num2str(k)]);
        set(gca, 'FontSize', fontsize)
        
        subplot(2, 1, 2);
        fe.plotWithSettings(mesh.nodes,"deformed",analysisSecondOrder.fromFEMVector( analysisSecondOrder.modes5(:,k+1) ),0.1,"elem nums",topOptSecondOrder.allx(:,k+1)>=0.5,"nodes",false);
        title(['Mode 5, vol_{fr}=' num2str(topOptSecondOrder.plVol(k+1)) ', Frq. =' num2str(topOptSecondOrder.plOmegas(5,k+1),4) ' [Hz]' ', iter' num2str(k+1)]);
        set(gca, 'FontSize', fontsize)
        saveas(gcf,['frame_correlation_mode_5_' num2str(k) '.pdf'])
        savefig(gcf,['frame_correlation_mode_5_' num2str(k) '.fig'])
    end
end

% figure, hold on
% p2=plot(topOptBuckling.plVol,topOptBuckling.plOmegas,'r','LineWidth', 3);
% set(gca, 'XDir', 'reverse');
% title('First 5 frequencies evolution');
% xlabel('Volume fracion [%]');
% ylabel('Frequency [Hz]');
% %xlim([37 57]);
% set(gca, 'FontSize', 24)
% 
% figure, hold on
% p2=plot(topOptBuckling.plVol,topOptBuckling.plLambda,'r','LineWidth', 3);
% set(gca, 'XDir', 'reverse');
% title('Critical force coefficient evolution with buckling');
% xlabel('Volume fracion [%]');
% ylabel('Critical force coefficient [%]');
% %xlim([37 57]);
% set(gca, 'FontSize', 24)

% figure, hold on
% p1=plot(topOptSecondOrder.plVol,topOptSecondOrder.plLambda,'b','LineWidth', 3);
% p2=plot(topOptBuckling.plVol,topOptBuckling.plLambda,'r','LineWidth', 3);
% legend([p1, p2], {'Without buckling', 'With buckling'});
% set(gca, 'XDir', 'reverse');
% title('Comparison of critical force evolution coefficient');
% xlabel('Volume fracion [%]');
% ylabel('Critical force coefficient [%]');
% %xlim([37 57]);
% set(gca, 'FontSize', 24)

% figure;
% tic
% topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, analysis, penal, 0.2, true);
% [objF, xopt]  = topOpt.solve();
% toc

%save('Cantilever2DBucklingDown80.mat');

