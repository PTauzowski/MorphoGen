function  plotArmTopologies(matFile,gtl)

load(matFile);

figure, hold on;
subplot(1,2,1);
topOptSecondOrder.plotCurrentFrame();
title('Without buckling constraints');

subplot(1,2,2);
topOptBuckling.plotCurrentFrame();
title('Including buckling constraints')
sgtitle(gtl);

end

