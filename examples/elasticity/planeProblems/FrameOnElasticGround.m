
clear;
close all;

Ecl  = 25.0E9; 
Ebm  = 25.0E9;
Egr  = 80.0E9;
nugr = 0.3;

b=0.25;
hbm=0.8;
hcl=0.7;
lspan=5;
hfloor=3;
nspan=3;
nfloor=3;


lgr=30;
hgr=15;
resgr=30;


model = FrameOnElasticGroundModel( Ecl, hcl, Ebm, hbm, b, lspan, hfloor, nspan, nfloor, Egr, nugr, lgr, hgr, resgr, ShapeFunctionL4 );
