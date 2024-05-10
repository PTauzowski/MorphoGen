
Vend=0.4;
x=(0:0.01:1)*Vend;
a=10;
b=1;
c=0.05;
y1=1./(1+exp(-a*(x/Vend-b)))-1./(1+exp(-a*(0-b))) + c;
y2=1-(1./(1+exp(-a*(x/Vend-b)))-1./(1+exp(-a*(0-b))));

%plot(x,y1);
plot(x,y2);
