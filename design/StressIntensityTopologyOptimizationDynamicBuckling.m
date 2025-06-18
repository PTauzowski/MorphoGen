classdef StressIntensityTopologyOptimizationDynamicBuckling < StressIntensityTopologyOptimization
    
    properties
         Vend, plLambda, plVol, plOmegas, lastStableFrame, bucklingForms, vibrationForms;
    end
    
    methods
        function obj = StressIntensityTopologyOptimizationDynamicBuckling(Rmin,linearElasticProblem,maxais,penal,Vend,is_const)
            obj=obj@StressIntensityTopologyOptimization(1,Rmin,linearElasticProblem,maxais,penal,is_const)
            obj.Vend=Vend;
            obj.plLambda=[];
            obj.plVol=[];
            obj.plOmegas=[];
        end
                           
        function of = computeObjectiveFunction(obj)
            of = sum( obj.x );
            obj.FobjValue=of;
        end

        function dc = computeInequalityConstraints(obj,x)
            dc = sum( x )/obj.V0 - obj.Vend;
        end

        function dc = computeEqualityConstraints(obj)
        end
        
        function printIterationInfo(obj)
            fprintf('%5i ',obj.iteration);
            fprintf('Vrel=%2.1f ',round(sum( obj.x )/obj.V0*1000)/10);
            %fprintf('lambda=%5.3g ', obj.FEAnalysis.lambda);
            fprintf('frq(1)=%5.3g ', obj.FEAnalysis.frequencies(1,end));
            fprintf('frq(2)=%5.3g ', obj.FEAnalysis.frequencies(2,end));
            fprintf('frq(3)=%5.3g ', obj.FEAnalysis.frequencies(3,end));
            fprintf('frq(4)=%5.3g ', obj.FEAnalysis.frequencies(4,end));
            fprintf('\n');
%            obj.plLambda = [ obj.plLambda abs( obj.FEAnalysis.lambda) ];
            obj.plOmegas = [ obj.plOmegas abs( obj.FEAnalysis.omegas(:,end) ) ];
            obj.plVol = [ obj.plVol round(sum( obj.x )/obj.V0*1000)/10 ];
        end

        function [plCor] = getModesSelfCorrelation(obj, mode)
           plCor=[];
           for k=2:size(mode,2)
               plCor = [plCor abs(mode(:,k-1)'*mode(:,k))/norm(mode(:,k-1))/norm(mode(:,k))];
           end
        end

        function corrIdx = getModeSwitchMatrix(obj)
           nmodes=size(obj.FEAnalysis.omegas,1);
           corrIdx=zeros(nmodes,obj.iteration-1);
           coors=zeros(1,nmodes);
           for l=1:obj.iteration-1
                for i=1:nmodes
                    for j=1:nmodes
                        coors(j)= abs(obj.FEAnalysis.modes(:,i,1)'*obj.FEAnalysis.modes(:,j,l))/norm(obj.FEAnalysis.modes(:,i,1))/norm(obj.FEAnalysis.modes(:,j,l));
                    end
                    [~, corrIdx(i,l)]=max(coors);
               end
           end
        end

        function plot_frequencies_switched(obj, basename, load_mode, nmodes, corrIdx)
            niter=size(obj.plOmegas,2);
            for l=1:nmodes
                figure, hold on
                o=zeros(1,niter);
                for i=1:niter
                    o(i) = obj.plOmegas(corrIdx(l,i),i)/2/pi;
                end
                p2=plot(obj.plVol,o','LineWidth', 3);
                set(gca, 'XDir', 'reverse');
                title(['Frequency with switch' newline 'Load mode '  num2str(load_mode)  ', Eigen mode '  num2str(l)]); 
                xlabel('Volume fracion [%]');
                ylabel('Frequency with switch [Hz]');
                ylim([0 inf])
                set(gca, 'FontSize', 18)
                saveas(gcf,[basename '_switched_frequency_loadmode_' num2str(load_mode) '_mode_'  num2str(l)  '.png']);
                savefig(gcf,[basename '_switched_frequency_loadmode_' num2str(load_mode) '_mode_'  num2str(l)  '.fig'])
            end
        end

        function plot_frequencies(obj, basename, load_mode, nmodes)
            for l=1:nmodes
                figure, hold on
                p2=plot(obj.plVol,obj.plOmegas(l,:)/2/pi','LineWidth', 3);
                set(gca, 'XDir', 'reverse');
                title(['Load mode '  num2str(load_mode)  ', Eigen mode '  num2str(l)  ' evolution']); 
                xlabel('Volume fracion [%]');
                ylabel('Frequency [Hz]');
                ylim([0 inf])
                set(gca, 'FontSize', 18)
                saveas(gcf,[basename '_frequency_loadmode_' num2str(load_mode) '_mode_'  num2str(l)  '.png']);
                savefig(gcf,[basename '_frequency_loadmode_' num2str(load_mode) '_mode_'  num2str(l)  '.fig'])
            end
        end

        function plot_forms(obj, basename, nmodes, frame, load_mode, scales)
            fontsize=18;
            fe=obj.FEAnalysis.felems{1};
            mesh=obj.FEAnalysis.mesh;
            for i=1:nmodes
                figure, hold on;
                %subplot(nmodes, 1, i);
                fe.plotWithSettings(mesh.nodes,"deformed",obj.FEAnalysis.fromFEMVector( obj.FEAnalysis.modes(:,i,frame) ),scales(i),"elem nums",obj.allx(:,frame)>=0.5,"edge color", "k");
                title(['Load mode ' num2str(load_mode) ', Eigen mode ' num2str(i) ', vol_{fr}=' num2str(obj.plVol(frame)) ', Frq. =' num2str(obj.plOmegas(i,frame)/2/pi,4) ' [Hz]'  ', iter:' num2str(frame)]);
                set(gca, 'FontSize', fontsize)
                exportgraphics(gcf,[basename '_loadmode_' num2str(load_mode) '_mode_' num2str(i)  '_iters_' num2str(frame) '.pdf'],'ContentType','vector')
                saveas(gcf,[basename '_loadmode_' num2str(load_mode) '_mode_' num2str(i)  '_iters_' num2str(frame) '.png']);
                savefig(gcf,[basename '_loadmode_' num2str(load_mode) '_mode_' num2str(i)  '_iters_' num2str(frame) '.fig'])
            end
           
        end

        function plot_subsequent_forms(obj, basename, mode1, mode2, frame)
            fontsize=16;

            figure;
            subplot(2, 1, 1);
            fe.plotWithSettings(mesh.nodes,"deformed",obj.FEAnalysis.fromFEMVector( obj.FEAnalysis.modes(:,mode1,frame) ),0.2,"elem nums",obj.allx(:,frame)>=0.5);
            title(['Mode ' num2str(mode1) ', vol_{fr}=' num2str(topOptSecondOrder.plVol(frame)) ', Frq. =' num2str(topOptSecondOrder.plOmegas(mode1,frame)/2/pi,4) ' [Hz]'  ', iter:' num2str(frame)]);
            set(gca, 'FontSize', fontsize)
            
            subplot(2, 1, 2);
            fe.plotWithSettings(mesh.nodes,"deformed",obj.FEAnalysis.fromFEMVector( obj.FEAnalysis.modes(:,mode2,frame) ),0.2,"elem nums",obj.allx(:,frame)>=0.5);
            title(['Mode ' num2str(mode2) ', vol_{fr}=' num2str(topOptSecondOrder.plVol(frame)) ', Frq. =' num2str(topOptSecondOrder.plOmegas(mode2,frame)/2/pi,4) ' [Hz]' ', iter' num2str(frame)]);
            set(gca, 'FontSize', fontsize)
            saveas(gcf,[ basename '_corrframes_' num2str(frame) '.png'])
            savefig(gcf,[ basename '_corrframes_' num2str(frame) '.fig'])

        end

        function plot_correlation_map(obj, basename, mode, frame)
            figure;
            U1 = normalize(obj.FEAnalysis.modes(:,mode1,frame-1), 1); 
            U2 = normalize(obj.FEAnalysis.modes(:,mode1,frame), 1);
            C = abs(U1' * U2);
            imagesc(C);
            colormap(flipud(hot));
            colorbar;
            axis equal tight;
            xlabel(['Mode number' num2str(mode)]);
            ylabel(['Frame ' num2str(frame)]);
            title('Correlation matrix of eigenmode 1');
            saveas(gcf,[ basename '_corr_map_' num2str(frame) '.png'])
            savefig(gcf,[ basename '_corrf_map_' num2str(frame) '.fig'])
        end

        function plot_correlation_curve(obj, basename, load_mode, nforms)
            figure;
            tol=1.0E3;
            for k=1:nforms
                kstr=num2str(k);
                %plCorr = round(obj.getModesCorrelation(squeeze(obj.FEAnalysis.modes(:,k,:)))/tol)*tol;
                plCorr = obj.getModesSelfCorrelation(squeeze(obj.FEAnalysis.modes(:,k,:)));
                figure, hold on
                p2=plot(obj.plVol(1:end-1), plCorr(1,1:end)','LineWidth', 2);
                ylim([0 1])
                set(gca, 'XDir', 'reverse');
                title(['Self correlations of mode ' kstr ' for load form' num2str(load_mode)]);
                xlabel('Volume fracion [%]');
                ylabel('Correlation');
                %xlim([37 57]);
                set(gca, 'FontSize', 18)
                saveas(gcf,[basename '_loadmode_' num2str(load_mode) '_correlation_' kstr '.png'])
                savefig(gcf,[basename '_loadmode_' num2str(load_mode) '_correlation_' kstr '.fig'])
            end
        end

        function plot_uncorrelated_frames(obj, basename)
            for mode=1:5
                kstr=num2str(mode);
                plCorr = obj.getModesSelfCorrelation(obj.FEAnalysis.modes(:,mode,:));
                for k=1:size(plCorr,2)
                    if plCorr(k)<0.3
                        figure;
                
                        subplot(2, 1, 1);
                        fe.plotWithSettings(mesh.nodes,"deformed",analysisSecondOrder.fromFEMVector( obj.FEAnalysis.modes(:,mode,k) ),0.1,"elem nums",topOptSecondOrder.allx(:,k)>=0.5,"nodes",false);
                        title(['Mode ' kstr ', vol_{fr}=' num2str(topOptSecondOrder.plVol(k)) ', Frq. =' num2str(topOptSecondOrder.plOmegas(1,k-1)/2/pi,4) ' [Hz]'], [ 'MAC=' num2str(pl_mode1_cor(k),3) ', iter:' num2str(k)]);
                        set(gca, 'FontSize', fontsize)
                        
                        subplot(2, 1, 2);
                        fe.plotWithSettings(mesh.nodes,"deformed",analysisSecondOrder.fromFEMVector( analysisSecondOrder.modes1(:,mode,k+1) ),0.1,"elem nums",topOptSecondOrder.allx(:,k+1)>=0.5,"nodes",false);
                        title(['Mode ' kstr ', vol_{fr}=' num2str(topOptSecondOrder.plVol(k+1)) ', Frq. =' num2str(topOptSecondOrder.plOmegas(1,k+1)/2/pi,4) ' [Hz]' ', iter' num2str(k+1)]);
                        set(gca, 'FontSize', fontsize)
                
                        saveas(gcf,[basename '_frame_correlation_mode_1_' num2str(k) '.png'])
                        savefig(gcf,[basename '_frame_correlation_mode_1_' num2str(k) '.fig'])
                    end
                end
            end
        end

        function prepareOutputImages(obj,nmodes)
            obj.plot_frequencies(nmodes);
            obj.plot_forms(5, 1);
            obj.plot_forms(5, size(obj.allx,2));
            obj.plot_correlation_curve();
            obj.plot_uncorrelated_frames();
        end
        
    end
end

