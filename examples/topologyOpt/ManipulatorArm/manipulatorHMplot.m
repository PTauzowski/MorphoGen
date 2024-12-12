function  manipulatorHMplot(model,endPoints)

figure;
model.fe.plot(model.mesh.nodes);
model.analysis.plotCurrentLoad();
line(endPoints(:,1),endPoints(:,2),endPoints(:,3),Marker=".",Color='r',LineStyle='none');
title('Model for minimal Huber-Mises');


x=ones(model.analysis.getTotalElemsNumber(),1);
model.analysis.solveWeighted(x);
model.analysis.computeElementResults();
model.analysis.plotMaps(["sHM"],0.1);
model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);

end

