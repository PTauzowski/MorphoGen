
Vend=0.4;
x=(0:0.01:1)*Vend;
a=20;
b=0;
tc1=0.05;
tc2=0.005;
y1=abs(tc2-tc1)*(2./(1+exp(-(a*(x-b))))-1)+min(tc1,tc2);
y2=abs(tc2-tc1)*(1-(2./(1+exp(-(a*(x-b))))-1))+min(tc1,tc2);

a=20;
b=0.2;
y3=abs(tc2-tc1)*(1./(1+exp(-(a*(x-b)))))+min(tc1,tc2);
y4=abs(tc2-tc1)*(1-(1./(1+exp(-(a*(x-b))))))+min(tc1,tc2);
y5=(1 - x/Vend) * 0.005 + (x/Vend) * 0.05;
y6=(1 - x/Vend) * 0.05 + (x/Vend) * 0.005;

plot(x,y5,x,y6,x,y1,x,y2,x,y3,x,y4,'LineWidth', 2);
hLegend = legend('linear ascending','linear descending','convex ascending', 'concave descending', 'sigmoid ascending', 'sigmoid descending','Location', 'east','FontSize', 20);

hTitle = title('Cut threshold functions','FontSize', 22);
xlabel('Volume Fraction','FontSize', 22);
ylabel('Cut threshold value','FontSize', 22);


% Set the font size of the legend
set(hLegend, 'FontSize', 22);
set(gca, 'FontSize', 22);

