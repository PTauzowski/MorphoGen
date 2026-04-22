function plotStep5(fileName, scheme_title, R, model, topOpt, forces)
    figure, hold on;
    hold on, axis off; 
    daspect([1 1 1]);
    xlabel("x");
    ylabel("y");
    zlabel("z");
    light('Position', [-1 -2 5], 'Style', 'local');
    light('Position', [1 1 5], 'Style', 'infinite');
    gca.FontSize = 18;

    topOpt.FEAnalysis.felems{1}.plot(topOpt.FEAnalysis.mesh.nodes);
    topOpt.FEAnalysis.felems{1}.face_alpha = 0.1;
    topOpt.FEAnalysis.felems{1}.edge_color = 'None';
    topOpt.FEAnalysis.felems{1}.plot([topOpt.FEAnalysis.mesh.nodes(:,1) topOpt.FEAnalysis.mesh.nodes(:,2) -topOpt.FEAnalysis.mesh.nodes(:,3)]);
    topOpt.FEAnalysis.felems{1}.face_alpha = 1.0;
    topOpt.FEAnalysis.plotSupport();

    a=0.4;
    coord = [...
    -a   -a  0;
     a   -a  0;
     a    a  0;
    -a    a  0];

    idx = [1 2 3 4]';

    xc = coord(:,1);
    yc = coord(:,2);
    zc = coord(:,3);
    patch(xc(idx), yc(idx), zc(idx), 'k', 'facealpha', 0.1);
    view(45,25);

    upward_facing_nodes=topOpt.FEAnalysis.mesh.findUpwardFacingNodes();
    nnodes=numel(upward_facing_nodes);

    load=zeros(nnodes,3);

    forces = forces/max(abs(forces))*0.2;

    N  = forces(1);
    Ty = forces(2);
    Tz = forces(3);

    Ms = forces(4);
    My = forces(5);
    Mz = forces(6);

    x=topOpt.FEAnalysis.mesh.nodes(upward_facing_nodes,:);

    xs=[0 0 0];
    
    load = load - [Tz 0 0];
    load = load + [0 Ty 0];
    load = load - [0 0 N];
    load = load - [x(:,1)*0 x(:,1)*0 (x(:,1)-xs(1))/R*My];
    load = load + [x(:,1)*0 x(:,1)*0 (x(:,2)/R)*Mz];
    load = load + Ms*[-x(:,2)./sqrt((x(:,1)-xs(1)).^2+x(:,2).^2) (x(:,1)-xs(1))./sqrt((x(:,1)-xs(1)).^2+x(:,2).^2) x(:,3)*0]; 

    X=topOpt.FEAnalysis.mesh.nodes(upward_facing_nodes,1)+load(:,1);
    Y=topOpt.FEAnalysis.mesh.nodes(upward_facing_nodes,2)+load(:,2);
    Z=topOpt.FEAnalysis.mesh.nodes(upward_facing_nodes,3)+load(:,3);

    % U=topOpt.FEAnalysis.mesh.nodes(upward_facing_nodes,1)-load(:,1);
    % V=topOpt.FEAnalysis.mesh.nodes(upward_facing_nodes,2)-load(:,2);
    % W=topOpt.FEAnalysis.mesh.nodes(upward_facing_nodes,3)-load(:,3);

    U=-load(:,1);
    V=-load(:,2);
    W=-load(:,3);

    quiver3(X,Y,Z,U,V,W,'LineWidth', 1,'AutoScale', 'off','MaxHeadSize',0.5,'Color','m');
    title(scheme_title);

    exportgraphics(gcf, fileName+".png", 'Resolution', 600); 
    savefig(gcf, fileName+".fig");
     
end