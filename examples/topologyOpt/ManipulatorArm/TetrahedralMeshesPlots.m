clear;
close all;
col=[0.8 0.8 0.8];

figure;
hold on;
daspect([1 1 1]);
load('TopologyTwoRingsEnvelope_tetrahedral_mesh.mat');
tetramesh(elems,nodes);
title('Tetrahedral two rings mesh of envelope topology using tetramesh');

figure;
hold on;
daspect([1 1 1]);
[~,ifaces,~] = unique( sort(elems,2),'rows' );
patch('Vertices', nodes, 'Faces', elems(ifaces,:),'FaceColor',col,'EdgeColor','k',"FaceAlpha",1);
title('Tetrahedral two rings mesh of envelope topology using patch');

figure;
hold on;
daspect([1 1 1]);
load('TopologyTwoRingsAverage_tetrahedral_mesh.mat');
tetramesh(elems,nodes);
title('Tetrahedral two rings verage topology using tetramesh');

figure;
hold on;
daspect([1 1 1]);
[~,ifaces,~] = unique( sort(elems,2),'rows' );
patch('Vertices', nodes, 'Faces', elems(ifaces,:),'FaceColor',col,'EdgeColor','k',"FaceAlpha",1);
title('Tetrahedral two rings mesh of average topology using patch');

figure;
hold on;
daspect([1 1 1]);
load('TopologyOneRingEnvelope_tetrahedral_mesh.mat');
tetramesh(elems,nodes);
title('Tetrahedral mesh of envelope topology using tetramesh');


figure;
hold on;
daspect([1 1 1]);
[~,ifaces,~] = unique( sort(elems,2),'rows' );
patch('Vertices', nodes, 'Faces', elems(ifaces,:),'FaceColor',col,'EdgeColor','k',"FaceAlpha",1);
title('Tetrahedral one ring mesh of envelope topology using patch');


figure;
hold on;
daspect([1 1 1]);
load('TopologyOneRingAverage_tetrahedral_mesh.mat');
tetramesh(elems,nodes);
title('Tetrahedral one ring mesh of average topology using tetramesh');


figure;
hold on;
daspect([1 1 1]);
[~,ifaces,~] = unique( sort(elems,2),'rows' );
patch('Vertices', nodes, 'Faces', elems(ifaces,:),'FaceColor',col,'EdgeColor','k',"FaceAlpha",1);
title('Tetrahedral one ring mesh of average topology using patch');