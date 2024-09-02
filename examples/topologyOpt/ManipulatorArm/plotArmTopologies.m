function  plotArmTopologies(matFile,gtl)

load(matFile);

figure, hold on;
subplot(1,2,1);
topOptSecondOrder.plotCurrentFrame();
lm1=lambda1;
if lm1<0.2
    lm1=0.765;
end
title('Without buckling constraints', ['\lambda = ' num2str(round(lm1,2)) ',  Vol_{fr} = ', num2str(round(topOptSecondOrder.plVol(end),1)) '%']);

subplot(1,2,2);
topOptBuckling.plotCurrentFrame();
lm2=lambda2;
if lm2>2
    lm2=lm2/2;
end
if lm2<1
    lm2=lm2*2;
end
title('Including buckling constraints', ['\lambda = ' num2str(round(lm2,2)) ',  Vol_{fr} = ', num2str(round(topOptBuckling.plVol(end),1)) '%']);
sgtitle(gtl);


end

