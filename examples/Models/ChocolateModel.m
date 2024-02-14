classdef ChocolateModel < ModelLinear

    properties
        ganTh, alGanTh, intTh, cx, cy, notchWidth;
    end
    
    methods
        function obj = ChocolateModel( ganTh, alGanTh, notchWidth, relNotchDepth, relRoutndNotchDepth, E, nu, alphaT, dT)
            obj.generateMesh( ganTh, alGanTh, notchWidth, relNotchDepth, relRoutndNotchDepth );
            obj.fe = SolidElasticElem( ShapeFunctionL27, obj.mesh.elems );
            obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, false );

            allganElemsSelector = Selector( @(x)( x(:,3) > ganTh ) );
            allganTopElemsSelector = Selector( @(x)( x(:,3) > ganTh+alGanTh*0.6 ) );
            %fixedFaceSelector = Selector( @(x)( abs(x(:,3) - l)<0.001 ) );
            %loadedFaceSelector = Selector( @(x)( abs(x(:,1) - l)<0.001 ) );

           
            meshMax=max(obj.mesh.nodes);
            obj.analysis.fixClosestNode([0 0 0], ["ux" "uy" "uz"], [0 0 0] );
            obj.analysis.fixClosestNode([meshMax(1) 0 0], ["uz"], 0);
            obj.analysis.fixClosestNode([meshMax(1) meshMax(2) ganTh], ["ux" "uz"], [0 0] );

            obj.analysis.loadElementsThermal(allganElemsSelector,alphaT*dT);
            obj.analysis.loadElementsThermal(allganTopElemsSelector,alphaT*0.4*dT);
            %obj.analysis.loadClosestNode([meshMax(1)/2 meshMax(2)/2 meshMax(2)], ["uz"], -1);
           
            obj.fe.props.h=1;
            material = SolidMaterial('mat1');
            material.setElasticIzo(E, nu);
            material.setElasticIzoGrad();
            obj.fe.setMaterial(material);
            
            obj.x=ones(1,obj.analysis.getTotalElemsNumber());
            obj.result_number=13;            
        end

        function generateMesh( obj, ganTh, alGanTh, notchWidth, relNotchDepth, relRoutndNotchDepth )
            obj.ganTh=ganTh;
            obj.alGanTh=alGanTh;
            obj.notchWidth=notchWidth;
            
            % GaN layer thickness
            %ganTh=5;
            
            % AlGanLayer thickness
            %alGanTh = 10;
            
            % interface thickness
            obj.intTh = 0.025;
            
            % width of the notch 
            %notchWidth=4; 
            
            %depth of the notch
            notchDepth=relNotchDepth*alGanTh; %0.875
            
            %depth of the part ow the notch which is rounded
            roundNotchDepth=relRoutndNotchDepth*notchDepth; %0.3
            
            %depth of the part ow the notch which is straight  
            straightNotchDepth=notchDepth-roundNotchDepth;
            
            % number of tiles in the x direction
            xtiles=4;
            
            % number of tiles in the y direction
            ytiles=4;
            
            % x - width of the tile
            cx=30-notchWidth;
            %cx=20-notchWidth;
            obj.cx=cx;
            
            % y - width of the tile
            cy=20-notchWidth; 
            obj.cy=cy;
            
            % FE x - division of the tile
            ncx=8;
            
            % FE y - division of the tile
            ncy=6;
            
            % depth FE division of the GaN layer
            ngan=2;
            
            % FE depth division of the rouned part of the notch
            nround=1;
            
            % FE depth division of the straight part of the notch
            nstr=2;
            
            % FE  width division of the half notch (second half is symetrical if inside)
            nnotch=2;
            
            % value of the chemistry imposed
            chemistry=0.1;
            
            % generated FEAP input file name 
            FEAPfilename = "chocolate4gan2.i";
            
            % One tile model generation
            ShapeFn = ShapeFunctionL8;
            ShapeFn27 = ShapeFunctionL27;
            mesh = Mesh();
            mesh2 = Mesh();
            % Xg1_8=[0 0 0;  cx+notchWidth 0 0; notchWidth/2 notchWidth/2 0; cx+notchWidth/2 notchWidth/2 0; ...
            %        0 0 notchDepth/4;  cx+notchWidth 0 notchDepth/4; notchWidth/2 notchWidth/2 notchDepth/2; cx+notchWidth/2 notchWidth/2 notchDepth/2]; 
            % 
            % Xg2_8=[ cx+notchWidth 0 0;  cx+notchWidth cy+notchWidth 0;  cx+notchWidth/2 notchWidth/2 0;   cx+notchWidth/2 cy+notchWidth/2 0; ...
            %         cx+notchWidth 0 notchDepth/4;  cx+notchWidth cy+notchWidth notchDepth/4;  cx+notchWidth/2 notchWidth/2 notchDepth/2;   cx+notchWidth/2 cy+notchWidth/2 notchDepth/2];
            
            Xg1_8=[0 0 0;  cx+notchWidth 0 0; notchWidth/2 notchWidth/2 0; cx+notchWidth/2 notchWidth/2 0; ...
                    0 0 notchDepth/2;  cx+notchWidth 0 notchDepth/2; notchWidth/2 notchWidth/2 notchDepth; cx+notchWidth/2 notchWidth/2 notchDepth]; 
             
            Xg2_8=[ cx+notchWidth 0 0;  cx+notchWidth cy+notchWidth 0;  cx+notchWidth/2 notchWidth/2 0;   cx+notchWidth/2 cy+notchWidth/2 0; ...
                     cx+notchWidth 0 notchDepth/2;  cx+notchWidth cy+notchWidth notchDepth/2;  cx+notchWidth/2 notchWidth/2 notchDepth; cx+notchWidth/2 cy+notchWidth/2 notchDepth/2];
            
            gh=alGanTh-notchDepth;
            gs=gh+roundNotchDepth/8;
            Xg1_27=[0 0 0;  (cx+notchWidth)/2 0 0; cx+notchWidth 0 0; ...
                    notchWidth/4 notchWidth/4 0; (cx+notchWidth)/2 notchWidth/4 0;  cx+3*notchWidth/4 notchWidth/4 0; ...
                    notchWidth/2 notchWidth/2 0; (cx+notchWidth)/2 notchWidth/2 0;  cx+notchWidth/2 notchWidth/2 0; ...
                    0 0 gh/2;  (cx+notchWidth)/2 0 gh/2; cx+notchWidth 0 gh/2;...
                    notchWidth/4 notchWidth/4 gs/2; (cx+notchWidth)/2 notchWidth/4 gs/2;  cx+3*notchWidth/4 notchWidth/4 gs/2; ...
                    notchWidth/2 notchWidth/2 roundNotchDepth/2; (cx+notchWidth)/2 notchWidth/2 roundNotchDepth/2;  cx+notchWidth/2 notchWidth/2 roundNotchDepth/2; ...
                    0 0 gh;  (cx+notchWidth)/2 0 gh; cx+notchWidth 0 gh;...
                    notchWidth/4 notchWidth/4 gs; (cx+notchWidth)/2 notchWidth/4 gs;  cx+3*notchWidth/4 notchWidth/4 gs; ...
                    notchWidth/2 notchWidth/2 roundNotchDepth; (cx+notchWidth)/2 notchWidth/2 roundNotchDepth;  cx+notchWidth/2 notchWidth/2 roundNotchDepth]; 
            
            Xg2_27=[cx+notchWidth   0              0;     cx+notchWidth      (cy+notchWidth)/2  0;  cx+notchWidth      cy+notchWidth  0; ...
                    cx+3*notchWidth/4 notchWidth/4 0;     cx+3*notchWidth/4  (cy+notchWidth)/2  0;  cx+3*notchWidth/4   cy+3*notchWidth/4    0; ...
                    cx+notchWidth/2 notchWidth/2   0;     cx+notchWidth/2    (cy+notchWidth)/2    0;  cx+notchWidth/2   cy+notchWidth/2    0; ...
                    cx+notchWidth   0              gh/2;  cx+notchWidth      (cy+notchWidth)/2  gh/2;  cx+notchWidth      cy+notchWidth  gh/2; ...
                    cx+3*notchWidth/4 notchWidth/4 gs/2;  cx+3*notchWidth/4  (cy+notchWidth)/2  gs/2;  cx+3*notchWidth/4   cy+3*notchWidth/4    gs/2; ...
                    cx+notchWidth/2 notchWidth/2   roundNotchDepth/2;  cx+notchWidth/2    (cy+notchWidth)/2    roundNotchDepth/2;  cx+notchWidth/2   cy+notchWidth/2    roundNotchDepth/2; ...
                    cx+notchWidth   0              gh;  cx+notchWidth      (cy+notchWidth)/2  gh;  cx+notchWidth      cy+notchWidth  gh; ...
                    cx+3*notchWidth/4 notchWidth/4 gs;  cx+3*notchWidth/4  (cy+notchWidth)/2  gs;  cx+3*notchWidth/4   cy+3*notchWidth/4    gs; ...
                    cx+notchWidth/2 notchWidth/2   roundNotchDepth;  cx+notchWidth/2    (cy+notchWidth)/2    roundNotchDepth;  cx+notchWidth/2   cy+notchWidth/2    roundNotchDepth]; 
            
            
            mesh.addShapedMesh3D( ShapeFn27, Xg1_27, [ncx,nnotch,nround], ShapeFn27.localNodes );
            mesh.duplicateTransformedMeshDeg3D( [(cx+notchWidth)/2  (cy+notchWidth)/2 ], 180, [0 0 0] );
            mesh2.addShapedMesh3D( ShapeFn27, Xg2_27, [ncy,nnotch,nround], ShapeFn27.localNodes );
            mesh2.duplicateTransformedMeshDeg3D( [(cx+notchWidth)/2  (cy+notchWidth)/2 ], 180, [0 0 0] );
            mesh.mergeMesh(mesh2);
            mesh.addRectMesh3D(notchWidth/2,notchWidth/2,0,cx,cy,roundNotchDepth,ncx,ncy,nround,ShapeFn27.localNodes);
            mesh.addRectMesh3D(notchWidth/2,notchWidth/2,roundNotchDepth,cx,cy,straightNotchDepth,ncx,ncy,nstr,ShapeFn27.localNodes);
            mesh.nodes=mesh.nodes+[0 0 obj.intTh];
            mesh.addRectMesh3D(notchWidth/2,notchWidth/2,0, cx,cy,obj.intTh,ncx,ncy,1,ShapeFn27.localNodes);
            
            Xg1a_27=[0 0 0;  (cx+notchWidth)/2 0 0; cx+notchWidth 0 0; ...
                     notchWidth/4 notchWidth/4 0; (cx+notchWidth)/2 notchWidth/4 0;  cx+3*notchWidth/4 notchWidth/4 0; ...
                     notchWidth/2 notchWidth/2 0; (cx+notchWidth)/2 notchWidth/2 0;  cx+notchWidth/2 notchWidth/2 0; ...
            
                     0 0 obj.intTh/2;  (cx+notchWidth)/2 0 obj.intTh/2; cx+notchWidth 0 obj.intTh/2;...
                     notchWidth/4 notchWidth/4 obj.intTh/2; (cx+notchWidth)/2 notchWidth/4 obj.intTh/2;  cx+3*notchWidth/4 notchWidth/4 obj.intTh/2; ...
                     notchWidth/2 notchWidth/2 obj.intTh/2; (cx+notchWidth)/2 notchWidth/2 obj.intTh/2;  cx+notchWidth/2 notchWidth/2 obj.intTh/2; ...
            
                     0 0 obj.intTh;  (cx+notchWidth)/2 0 obj.intTh; cx+notchWidth 0 obj.intTh;...
                     notchWidth/4 notchWidth/4 obj.intTh; (cx+notchWidth)/2 notchWidth/4 obj.intTh;  cx+3*notchWidth/4 notchWidth/4 obj.intTh; ...
                     notchWidth/2 notchWidth/2 obj.intTh; (cx+notchWidth)/2 notchWidth/2 obj.intTh;  cx+notchWidth/2 notchWidth/2 obj.intTh]; 
            
            Xg2a_27=[cx+notchWidth   0              0;     cx+notchWidth      (cy+notchWidth)/2  0;  cx+notchWidth      cy+notchWidth  0; ...
                    cx+3*notchWidth/4 notchWidth/4 0;     cx+3*notchWidth/4  (cy+notchWidth)/2  0;  cx+3*notchWidth/4   cy+3*notchWidth/4    0; ...
                    cx+notchWidth/2 notchWidth/2   0;     cx+notchWidth/2    (cy+notchWidth)/2    0;  cx+notchWidth/2   cy+notchWidth/2    0; ...
            
                    cx+notchWidth   0              obj.intTh/2;  cx+notchWidth      (cy+notchWidth)/2  obj.intTh/2;  cx+notchWidth      cy+notchWidth  obj.intTh/2; ...
                    cx+3*notchWidth/4 notchWidth/4 obj.intTh/2;  cx+3*notchWidth/4  (cy+notchWidth)/2  obj.intTh/2;  cx+3*notchWidth/4   cy+3*notchWidth/4    obj.intTh/2; ...
                    cx+notchWidth/2 notchWidth/2   obj.intTh/2;  cx+notchWidth/2    (cy+notchWidth)/2  obj.intTh/2;  cx+notchWidth/2   cy+notchWidth/2    obj.intTh/2; ...
             
                    cx+notchWidth   0              obj.intTh;  cx+notchWidth      (cy+notchWidth)/2  obj.intTh;  cx+notchWidth      cy+notchWidth  obj.intTh; ...
                    cx+3*notchWidth/4 notchWidth/4 obj.intTh;  cx+3*notchWidth/4  (cy+notchWidth)/2  obj.intTh;  cx+3*notchWidth/4   cy+3*notchWidth/4    obj.intTh; ...
                    cx+notchWidth/2 notchWidth/2   obj.intTh;  cx+notchWidth/2    (cy+notchWidth)/2  obj.intTh;  cx+notchWidth/2   cy+notchWidth/2    obj.intTh]; 
            
            mesh2=Mesh();
            mesh2.addShapedMesh3D( ShapeFn27, Xg1a_27, [ncx,nnotch,1], ShapeFn27.localNodes );
            mesh2.addShapedMesh3D( ShapeFn27, Xg2a_27, [ncy,nnotch,1], ShapeFn27.localNodes );
            mesh2.duplicateTransformedMeshDeg3D( [(cx+notchWidth)/2  (cy+notchWidth)/2 ], 180, [0 0 0] );
            mesh.mergeMesh(mesh2);
            
            
            mesh.nodes=mesh.nodes+[0 0 ganTh];
            mesh.addRectMesh3D(notchWidth/2,notchWidth/2,0, cx,cy,ganTh,ncx,ncy,ngan,ShapeFn27.localNodes);
            
            
            Xg1a_27=[0 0 0;  (cx+notchWidth)/2 0 0; cx+notchWidth 0 0; ...
                     notchWidth/4 notchWidth/4 0; (cx+notchWidth)/2 notchWidth/4 0;  cx+3*notchWidth/4 notchWidth/4 0; ...
                     notchWidth/2 notchWidth/2 0; (cx+notchWidth)/2 notchWidth/2 0;  cx+notchWidth/2 notchWidth/2 0; ...
            
                     0 0 ganTh/2;  (cx+notchWidth)/2 0 ganTh/2; cx+notchWidth 0 ganTh/2;...
                     notchWidth/4 notchWidth/4 ganTh/2; (cx+notchWidth)/2 notchWidth/4 ganTh/2;  cx+3*notchWidth/4 notchWidth/4 ganTh/2; ...
                     notchWidth/2 notchWidth/2 ganTh/2; (cx+notchWidth)/2 notchWidth/2 ganTh/2;  cx+notchWidth/2 notchWidth/2 ganTh/2; ...
            
                     0 0 ganTh;  (cx+notchWidth)/2 0 ganTh; cx+notchWidth 0 ganTh;...
                     notchWidth/4 notchWidth/4 ganTh; (cx+notchWidth)/2 notchWidth/4 ganTh;  cx+3*notchWidth/4 notchWidth/4 ganTh; ...
                     notchWidth/2 notchWidth/2 ganTh; (cx+notchWidth)/2 notchWidth/2 ganTh;  cx+notchWidth/2 notchWidth/2 ganTh]; 
            
            Xg2a_27=[cx+notchWidth   0              0;     cx+notchWidth      (cy+notchWidth)/2  0;  cx+notchWidth      cy+notchWidth  0; ...
                    cx+3*notchWidth/4 notchWidth/4 0;     cx+3*notchWidth/4  (cy+notchWidth)/2  0;  cx+3*notchWidth/4   cy+3*notchWidth/4    0; ...
                    cx+notchWidth/2 notchWidth/2   0;     cx+notchWidth/2    (cy+notchWidth)/2    0;  cx+notchWidth/2   cy+notchWidth/2    0; ...
            
                    cx+notchWidth   0              ganTh/2;  cx+notchWidth      (cy+notchWidth)/2  ganTh/2;  cx+notchWidth      cy+notchWidth  ganTh/2; ...
                    cx+3*notchWidth/4 notchWidth/4 ganTh/2;  cx+3*notchWidth/4  (cy+notchWidth)/2  ganTh/2;  cx+3*notchWidth/4   cy+3*notchWidth/4    ganTh/2; ...
                    cx+notchWidth/2 notchWidth/2   ganTh/2;  cx+notchWidth/2    (cy+notchWidth)/2  ganTh/2;  cx+notchWidth/2   cy+notchWidth/2    ganTh/2; ...
             
                    cx+notchWidth   0              ganTh;  cx+notchWidth      (cy+notchWidth)/2  ganTh;  cx+notchWidth      cy+notchWidth  ganTh; ...
                    cx+3*notchWidth/4 notchWidth/4 ganTh;  cx+3*notchWidth/4  (cy+notchWidth)/2  ganTh;  cx+3*notchWidth/4   cy+3*notchWidth/4    ganTh; ...
                    cx+notchWidth/2 notchWidth/2   ganTh;  cx+notchWidth/2    (cy+notchWidth)/2    ganTh;  cx+notchWidth/2   cy+notchWidth/2    ganTh]; 
            
            mesh2=Mesh();
            mesh2.addShapedMesh3D( ShapeFn27, Xg1a_27, [ncx,nnotch,ngan], ShapeFn27.localNodes );
            mesh2.addShapedMesh3D( ShapeFn27, Xg2a_27, [ncy,nnotch,ngan], ShapeFn27.localNodes );
            mesh2.duplicateTransformedMeshDeg3D( [(cx+notchWidth)/2  (cy+notchWidth)/2 ], 180, [0 0 0] );
            mesh.mergeMesh(mesh2);
            
            % copy tiles in x direction
            mesh.array(1,xtiles-1);
            
            % copy tiles in y direction
            mesh.array(2,ytiles-1);
            obj.mesh=mesh;
        end

        function s = getNodalStress(obj,result_node,result_number)
            sHM=obj.fe.results.nodal(result_node,result_number);
        end

        function [stressObj, sx1, sy1, sx2, sy2]=computeStressObjective(obj)
            x1=[(obj.cx+obj.notchWidth)/2 (obj.cy+obj.notchWidth)/2 0];
            x2=[obj.cx+obj.notchWidth obj.cy+obj.notchWidth 0];
            n1 = obj.mesh.findClosestNode(x1);
            n2 = obj.mesh.findClosestNode(x2);
            sx1 = obj.fe.results.nodal.all(n1,6);
            sy1 = obj.fe.results.nodal.all(n1,7);
            sx2 = obj.fe.results.nodal.all(n2,6);
            sy2 = obj.fe.results.nodal.all(n2,7);
            stressObj=sx1+sy1-(sx2+sy2);
            %stressObj=(sx2+sy2)-(sx1+sy1);
        end

        function o = computeObjectiveValue(obj,x)
            ganTh=h0*(1-x(1));
            alGanTh = h0*x(1);
            w=ganTh*x(2);
            r=h0*x(3);
            generateMesh( obj, ganTh, alGanTh, x(3), x(2) )
        end

        function FEAP_Export(obj,filename,mesh,ganTh,alGanTh,intTh,chemistry)
            myfile   = fopen ( filename, "w" );
            feapNum = [1 3 9 7 19 21 27 25 2 6 8 4 20 24 26 22 10 12 18 16  5 23 13 15 11 17 14];
            %feapNum = [1 9 2 12 21 10 4 11 3 17 25 18 23 27 24 20 26 19 5 13 6 16 22 14 8 15 7 ];
            fprintf( myfile, "feap * * chocolate \n  %d %d 0 3 4 27 \n\n", size(mesh.nodes,1), size(mesh.elems,1) );
            fprintf( myfile, "COORdinates\n");
            for k=1:size(mesh.nodes, 1)
                fprintf( myfile, "%d 0 %7.4E %7.4E  %7.4E\n",k,mesh.nodes(k,1),mesh.nodes(k,2),mesh.nodes(k,3));
            end
        
            fprintf( myfile, "\n\nELEMents\n" );
            for k=1:size(mesh.elems,1)
                felems = mesh.elems(k,feapNum);
                fprintf( myfile, "%d 0 1",k);
                for e=1:size(feapNum,2)
                    fprintf( myfile, " %d",felems(e) );  
                    if e==13
                        fprintf( myfile, "\n");
                    end
                end
                fprintf( myfile, "\n");
            end
            fprintf( myfile, "\n");
        
            zCoords = sort(unique(round(mesh.nodes(:,3)*10000)/10000));
        
            fprintf(myfile,"\n EDIS\n");
            fprintf(myfile,"  gap 0.001\n");
            nganTh=ganTh;
            for k=size(zCoords(:)):-1:1
                if ( zCoords(k) <= nganTh )
                    c=0;
                else
                    if abs( zCoords(k) - (nganTh + intTh/2)  ) < 1.0E-6 
                        c=chemistry/2;
                    elseif abs( zCoords(k) > ganTh+0.6*alGanTh  )      
                        c=1.4*chemistry;
                    else
                        c=chemistry;
                    end
                end
               fprintf(myfile,"  3  %5.3f  0  0  0  %1.2f\n", zCoords(k), c ); 
            end
                  
          meshMax=max(mesh.nodes);
          n1 = mesh.findClosestNode(  [0 0 0] );
          n2 = mesh.findClosestNode(  [meshMax(1) 0 0] );
          n3 = mesh.findClosestNode(  [meshMax(1) meshMax(2) ganTh] );
          xMax = zCoords(end);
          fprintf(myfile,"\n");  
           
          fprintf(myfile,"\n");
          fprintf(myfile,"BOUNdary\n");
          fprintf(myfile,"  1 1 0 0 0 -1\n");
          fprintf(myfile," %d 0 0 0 0  1\n",size(mesh.nodes,1));
          fprintf(myfile," %d 0 1 1 1  1\n",n1(1));
          fprintf(myfile," %d 0 0 0 1  1\n",n2(1));
          fprintf(myfile," %d 0 1 0 1  1\n",n3(1));
          fprintf(myfile,"\n");
        
          fprintf(myfile,"mate,1\n");
          fprintf(myfile,"user,14\n");
          fprintf(myfile,"3,2,2,1,30,0.,1.\n");
          fprintf(myfile,"2 -1 -1 0  0 0 0 1\n");
          fprintf(myfile,"3.189E-10 5.185E-10\n");
          fprintf(myfile,"div(sig) -2 1   1\n");
          fprintf(myfile,"0.\n");
          fprintf(myfile,"390.0d9,145.0d9,106.0d9,398.0d9,105.0d9\n");
          fprintf(myfile,"396.0d9,137.0d9,108.0d9,373.0d9,116.0d9\n");
          fprintf(myfile,"div(x_n) -1 1  -4\n");
          fprintf(myfile,"3.112E-10  4.982E-10\n");
          fprintf(myfile,"\n");
          fprintf(myfile,"end\n");
          fprintf(myfile,"\n");
          fprintf(myfile,"TIE\n");
          fprintf(myfile,"\n");
          fprintf(myfile,"BATCh\n");
          fprintf(myfile,"PROP\n");
          fprintf(myfile,"DT,,1\n");
          fprintf(myfile,"END\n");
          fprintf(myfile,"2 2\n");
          fprintf(myfile,"0 0 1 1\n");
          fprintf(myfile,"\n");
          fprintf(myfile,"batch\n");
          fprintf(myfile,"plot,pers,1\n");
          fprintf(myfile,"plot,hide\n");
          fprintf(myfile,"plot,fill\n");
          fprintf(myfile,"plot,defo\n");
          fprintf(myfile,"plot,mesh\n");
          fprintf(myfile,"plot,load\n");
          fprintf(myfile,"plot,axis\n");
          fprintf(myfile,"end\n");
          fprintf(myfile,"0\n");
          fprintf(myfile,"2000. 4000. 2000.\n");
          fprintf(myfile,"  0.  0. 2.\n");
          fprintf(myfile,"\n");
          fprintf(myfile,"batch\n");
          fprintf(myfile,"plot,mesh\n");
          fprintf(myfile,"plot,defo,1,1\n");
          fprintf(myfile,"end\n");
          fprintf(myfile,"\n");
          fprintf(myfile,"batch\n");
          fprintf(myfile,"opti\n");
          fprintf(myfile,"\n");
          fprintf(myfile,"LOOP,,1\n");
          fprintf(myfile,"TIME\n");
          fprintf(myfile,"LOOP,,99\n");
          fprintf(myfile,"utan,,1\n");
          fprintf(myfile,"plot,cont,4\n");
          fprintf(myfile,"plot,stre,1\n");
          fprintf(myfile,"plot,stre,2\n");
          fprintf(myfile,"plot,stre,3\n");
          fprintf(myfile,"plot,stre,4\n");
          fprintf(myfile,"plot,stre,5\n");
          fprintf(myfile,"!plot,stre,6\n");
          fprintf(myfile,"NEXT\n");
          fprintf(myfile,"save\n");
          fprintf(myfile,"disp,all\n");
          fprintf(myfile,"stre,all\n");
          fprintf(myfile,"stre,node,all\n");
          fprintf(myfile,"NEXT\n");
          fprintf(myfile,"\n");
          fprintf(myfile,"end\n");
          fprintf(myfile,"\n");
          fprintf(myfile,"inte\n");
          fprintf(myfile,"stop\n");
          fclose(myfile);
        end
    end
end

