clear;
close all;
col=[0.8 0.8 0.8];

edges = [1 2; 2 4; 4 1; 2 3; 3 1; 3 4 ]';
faces = [1 2 4; 2 3 4; 3 1 4; 2 1 3]';
fcontours = [1 2 4; 2 3 4; 3 1 4; 2 1 3]';

figure;
hold on;
daspect([1 1 1]);
load('TopologyTwoRingsEnvelope_tetrahedral_mesh.mat');
tetramesh(elems,nodes);
title('Tetrahedral two rings mesh of envelope topology using tetramesh');

figure;
hold on;
daspect([1 1 1]);
allfaces = reshape(elems(:,fcontours)',size(fcontours,1),size(fcontours,2)*size(elems,1))';
alledges = reshape(elems(:,edges)',size(edges,1),size(edges,2)*size(elems,1))';
[~,ifaces,~] = unique( sort(reshape(elems(:,fcontours)',size(fcontours,1),size(fcontours,2)*size(elems,1))',2),'rows' );
[uniqueEdges, ~, idx] = unique(sort(reshape(elems(:,edges)',size(edges,1),size(edges,2)*size(elems,1))',2), 'rows', 'stable');
counts = histcounts(idx, 1:(max(idx)+1));
edgesNoDuplicates = uniqueEdges(counts == 1, :);
patch('Vertices', nodes, 'Faces', allfaces(ifaces,:),'FaceColor',col,'EdgeColor','none',"FaceAlpha",1.0);
patch('Vertices', nodes, 'Faces', allfaces(ifaces,:),'FaceColor','none','EdgeColor','k');
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
allfaces = reshape(elems(:,fcontours)',size(fcontours,1),size(fcontours,2)*size(elems,1))';
alledges = reshape(elems(:,edges)',size(edges,1),size(edges,2)*size(elems,1))';
[~,ifaces,~] = unique( sort(reshape(elems(:,fcontours)',size(fcontours,1),size(fcontours,2)*size(elems,1))',2),'rows' );
[uniqueEdges, ~, idx] = unique(sort(reshape(elems(:,edges)',size(edges,1),size(edges,2)*size(elems,1))',2), 'rows', 'stable');
counts = histcounts(idx, 1:(max(idx)+1));
edgesNoDuplicates = uniqueEdges(counts == 1, :);
patch('Vertices', nodes, 'Faces', allfaces(ifaces,:),'FaceColor',col,'EdgeColor','none',"FaceAlpha",1.0);
patch('Vertices', nodes, 'Faces', allfaces(ifaces,:),'FaceColor','none','EdgeColor','k');
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
allfaces = reshape(elems(:,fcontours)',size(fcontours,1),size(fcontours,2)*size(elems,1))';
alledges = reshape(elems(:,edges)',size(edges,1),size(edges,2)*size(elems,1))';
[~,ifaces,~] = unique( sort(reshape(elems(:,fcontours)',size(fcontours,1),size(fcontours,2)*size(elems,1))',2),'rows' );
[uniqueEdges, ~, idx] = unique(sort(reshape(elems(:,edges)',size(edges,1),size(edges,2)*size(elems,1))',2), 'rows', 'stable');
counts = histcounts(idx, 1:(max(idx)+1));
edgesNoDuplicates = uniqueEdges(counts == 1, :);
patch('Vertices', nodes, 'Faces', allfaces(ifaces,:),'FaceColor',col,'EdgeColor','none',"FaceAlpha",1.0);
patch('Vertices', nodes, 'Faces', allfaces(ifaces,:),'FaceColor','none','EdgeColor','k');
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
allfaces = reshape(elems(:,fcontours)',size(fcontours,1),size(fcontours,2)*size(elems,1))';
alledges = reshape(elems(:,edges)',size(edges,1),size(edges,2)*size(elems,1))';
[~,ifaces,~] = unique( sort(reshape(elems(:,fcontours)',size(fcontours,1),size(fcontours,2)*size(elems,1))',2),'rows' );
[uniqueEdges, ~, idx] = unique(sort(reshape(elems(:,edges)',size(edges,1),size(edges,2)*size(elems,1))',2), 'rows', 'stable');
counts = histcounts(idx, 1:(max(idx)+1));
edgesNoDuplicates = uniqueEdges(counts == 1, :);
patch('Vertices', nodes, 'Faces', allfaces(ifaces,:),'FaceColor',col,'EdgeColor','none',"FaceAlpha",1.0);
patch('Vertices', nodes, 'Faces', allfaces(ifaces,:),'FaceColor','none','EdgeColor','k');
title('Tetrahedral one ring mesh of average topology using patch');