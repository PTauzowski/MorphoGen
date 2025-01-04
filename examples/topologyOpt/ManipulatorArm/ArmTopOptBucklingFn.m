function  [xopt, xopt_buckling, lambda1,lambda2] = ArmTopOptBucklingFn(name,R,r,h,alpha,Tx,Ty,N,Mx,My,Ms)

    th=R-r;

    
    % resCirc=100;
    % resTh=max(1,round(th/2/pi/R*resCirc));
    % resHeight=max(1,round(h/4/pi/R*resCirc));
    
    resTh=1;
    resCirc=round(2*pi*R/th*resTh);
    resHeight=round(h/th*resTh);
    
    % Filtering radius
    Rfilter = R*3*pi/resCirc;
    penal=3;
    cutTreshold = 0.005;
    
    ShapeFn = ShapeFunctionL8;
    mesh = Mesh();
    mesh.addManipulatorHalfSegment3D(r, R, h, alpha*pi/180, resTh, resCirc, resHeight, ShapeFn.localNodes);
    fe = SolidElasticElem( ShapeFn, mesh.elems );
    fe.plot(mesh.nodes);
    upward_facing_nodes=mesh.findUpwardFacingNodes();
    downward_facing_nodes=mesh.findDownwardFacingNodes();
    xs=(min(mesh.nodes)+max(mesh.nodes))/2;
    % scatter3(mesh.nodes(upward_facing_nodes,1), mesh.nodes(upward_facing_nodes,2), mesh.nodes(upward_facing_nodes,3), ...
    %          100, 'r', 'filled');
    
    fe.props.h=1;
    material = SolidMaterial('mat1');
    material.setElasticIzo(210.0E9, 0.3);
    material.setElasticIzoGrad();
    fe.setMaterial(material);
    
    analysis = LinearElasticityWeighted( fe, mesh, false );
    %problem = LinearElasticity( fe, mesh );
    % fixedEdgeSelector = Selector( @(x)( abs(x(:,3)) < 0.001 ) );
    % loadedFaceSelector = Selector( @(x)( abs(x(:,3)- Length) < 0.001 ) );
    constElemsSelector =  Selector( @(x)( (x(:,3) < 0.05 * h ) ) & (x(:,3) > 0.96 * h ) );
    
    
    un=false(size(mesh.nodes,1),1);
    dn=false(size(mesh.nodes,1),1);
    un( upward_facing_nodes ) = true;
    dn(downward_facing_nodes ) = true;
    loadedFaceSelector = Selector( un );
    fixedEdgeSelector = Selector( dn );
    upElems = mesh.findElems( loadedFaceSelector, false );
    downElems = mesh.findElems( fixedEdgeSelector, false );
    
    analysis.elementLoadSurfaceIntegral( "global", loadedFaceSelector, ["ux" "uy" "uz"], @(x)( x*0 + [Tx 0 0] ));
    analysis.elementLoadSurfaceIntegral( "global", loadedFaceSelector, ["ux" "uy" "uz"], @(x)( x*0 + [0 Ty 0] ));
    analysis.elementLoadSurfaceIntegral( "global", loadedFaceSelector, ["ux" "uy" "uz"], @(x)( x*0 + [0 0 N] ));
    analysis.elementLoadSurfaceIntegral( "global", loadedFaceSelector, ["ux" "uy" "uz"], @(x)( x*0 + [x(:,1)*0 x(:,1)*0 (x(:,1)-xs(1))*My] ));
    analysis.elementLoadSurfaceIntegral( "global", loadedFaceSelector, ["ux" "uy" "uz"], @(x)( x*0 + [x(:,1)*0 x(:,1)*0 (x(:,2)/R)*Mx] ));
    analysis.elementLoadSurfaceIntegral( "global", loadedFaceSelector, ["ux" "uy" "uz"], @(x)( x*0 + Ms*[-x(:,2)./sqrt((x(:,1)-xs(1)).^2+x(:,2).^2) (x(:,1)-xs(1))./sqrt((x(:,1)-xs(1)).^2+x(:,2).^2) x(:,3)*0] )); 
    
    analysis.fixNodes( fixedEdgeSelector, ["ux" "uy" "uz"] );
    %analysis.fixClosestNode( [0 0 0], ["ux" "uy" "uz"], [0 0 0]);
    %const_elems = [upElems downElems];
    const_elems = upElems;
    
    %mesh.transformNodesXY( @(x)( [ x(:,1) x(:,2) x(:,3)-0.3*x(:,1).*x(:,3)/Length ] )  );
    
    fe.plot(mesh.nodes);
    %problem.plotNodes();
    analysis.plotCurrentLoad();
    analysis.plotSupport();
    
    view(45, 45);
    
    analysis.printProblemInfo();
    
    nEigenForms=10;
    stability = LinearStability( analysis.felems, mesh);
    stability.Pnodal = analysis.Pnodal;
    stability.Pfem = analysis.Pfem;
    stability.supports = analysis.supports;
    stability.solve( nEigenForms);
    lambdas = diag(stability.lambdas)
    for k=1:min(2,nEigenForms)
        figure;
        %subplot(5, 2, k);
        stability.setForm(k);
        fe.plotSolidDeformed(mesh.nodes,stability.qnodal,0.2);
        axis on, xlabel('x-axis'), ylabel('y-axis'), view(3)
        lambda_str = sprintf('%.4g', lambdas(k));
        title(['Form:' num2str(k), ' \lambda=' lambda_str]);
    end
    
    analysisSecondOrder = SecondOrderElasticityWeighted( fe, mesh, 0.0, false );
    analysisSecondOrder.Pnodal=stability.Pnodal;
    analysisSecondOrder.Pfem=stability.Pfem;
    analysisSecondOrder.supports=stability.supports;
    
    analysisWithBuckling = SecondOrderElasticityWeighted( fe, mesh, 0.95, false );
    analysisWithBuckling.Pnodal=stability.Pnodal;
    analysisWithBuckling.Pfem=stability.Pfem;
    analysisWithBuckling.supports=stability.supports;
    
    % figure;
    % tic
    % topOptLinear = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penal, 0.2, false );
    % [objF, xopt]  = topOptLinear.solve();
    % toc
    
    % figure;
    % tic
    % topOptSecondOrder = StressIntensityTopologyOptimizationBuckling( Rfilter, analysisSecondOrder, cutTreshold, penal, 0.4, true );
    % topOptSecondOrder.setConstElems(const_elems);
    % [objF, xopt]  = topOptSecondOrder.solve();
    % toc
    
    figure;
    tic
    topOptBuckling = StressIntensityTopologyOptimizationBuckling( Rfilter, analysisWithBuckling, cutTreshold, penal, 0.4, true );
    topOptBuckling.setConstElems(const_elems);
    [objF, xopt]  = topOptBuckling.solve();
    toc
    savefig([name '_b.fig']);
    
    Pcr=abs(topOptBuckling.plLambda(end));
    Pfem=analysisWithBuckling.Pfem;
    stability.Pfem=Pfem*Pcr;
    stability.solveWeighted( xopt, nEigenForms );
    lambda2=abs(stability.lambdas(1));
    
    
   
    figure;
    tic
    analysisSecondOrder.Pfem=analysisSecondOrder.Pfem*Pcr;
    topOptSecondOrder = StressIntensityTopologyOptimizationBuckling( Rfilter, analysisSecondOrder, cutTreshold, penal, 0.4, true );
    topOptSecondOrder.setConstElems(const_elems);
    [objF, xopt_buckling]  = topOptSecondOrder.solve();
    lambda1=abs(topOptSecondOrder.plLambda(end));
    toc
    savefig([name '_l.fig']);
    
    % if abs(topOptBuckling.plLambda(end)) > 1
    %     vb=abs(topOptBuckling.plLambda(end));
    %     b=abs(topOptBuckling.plLambda(end));
    %     a=b;
    %     va=vb;
    %     while va>=1
    %         a=a*2;
    %         stability.Pfem=a*Pfem;  
    %         stability.solveWeighted( xopt, nEigenForms );
    %         va=abs(stability.lambdas(1));
    %     end
    % else
    %     va=abs(topOptBuckling.plLambda(end));
    %     a=abs(topOptBuckling.plLambda(end));
    %     b=a;
    %     vb=va;
    %     while vb<=1
    %         b=b/2;
    %         stability.Pfem=b*Pfem; 
    %         stability.solveWeighted( xopt, nEigenForms ); 
    %         vb=abs(stability.lambdas(1)) ;
    %     end
    % end
    % 
    % s=(a+b)/2;
    % stability.Pfem=s*Pfem;   
    % stability.solveWeighted( xopt, nEigenForms ); 
    % vs=abs(stability.lambdas(1));
    % while abs(b-a) > 0.0001
    %     if vs < 1
    %         a=s;
    %         va=vs;
    %     else
    %         b=s;
    %         vb=vs;
    %     end
    %     stability.Pfem=s*Pfem;    
    %     stability.solveWeighted( xopt, nEigenForms ); 
    %     vs=abs(stability.lambdas(1));
    %     s=(a+b)/2;
    %     fprintf('a=%5.4g, va=%5.4g, b=%5.4g, vb=%5.4g, s=%5.4g, vs=%5.4g\n', a,va,b,vb,s,vs);
    % end

     %fprintf('Topology: a=%s, Pcr=%5.4g, Lambda1=%5.4g, lambda2=%5.4g\n', name,topOptBuckling.plLambda(end),lambda1,lambda2);
     %save(['ArmTopopt_' name '.mat']);
end

