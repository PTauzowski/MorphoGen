classdef ChocolateModel < ModelLinear

    properties
        ganTh, alGanTh, intTh, tileXWidth, tileYWidth, notchWidth, nTempVars,xt,xn,zCoords,zCornersCoords;
        allganElemsSelector,allganTopElemsSelector,alphaT, dT, zTol;
    end
    
    methods
        function obj = ChocolateModel( ganTh, alGanTh, notchWidth, relNotchDepth, relRoutndNotchDepth, E, nu, alphaT, dT)
            obj.alphaT=alphaT;
            obj.dT=dT;
            obj.zTol=1000;
            obj.generateMesh( ganTh, alGanTh, notchWidth, relNotchDepth, relRoutndNotchDepth );
            obj.zCoords = sort(unique(round(obj.mesh.nodes(:,3)*obj.zTol)/obj.zTol),'descend');
            obj.computeZNodalCoords();
            obj.nTempVars=size(find((round((obj.zCoords-ganTh)*obj.zTol)/obj.zTol)>0),1)-3;
            obj.fe = SolidElasticElem( ShapeFunctionL27, obj.mesh.elems );
            obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, false );

            obj.allganElemsSelector = Selector( @(x)( x(:,3) > ganTh ) );
            obj.allganTopElemsSelector = Selector( @(x)( x(:,3) > ganTh+alGanTh*0.6 ) );
            %fixedFaceSelector = Selector( @(x)( abs(x(:,3) - l)<0.001 ) );
            %loadedFaceSelector = Selector( @(x)( abs(x(:,1) - l)<0.001 ) );

            meshMax=max(obj.mesh.nodes);
            obj.analysis.fixClosestNode([0 0 0], ["ux" "uy" "uz"], [0 0 0] );
            obj.analysis.fixClosestNode([meshMax(1) 0 0], ["uz"], 0);
            obj.analysis.fixClosestNode([meshMax(1) meshMax(2) ganTh], ["ux" "uz"], [0 0] );
           
            obj.fe.props.h=1;
            obj.fe.props.ndT=zeros(1,size(obj.mesh.nodes,1));
            material = SolidMaterial('mat1');
            material.setElasticIzo(E, nu);
            material.setElasticIzoGrad();
            obj.fe.setMaterial(material);

            %obj.analysis.loadClosestNode([meshMax(1)/2 meshMax(2)/2 meshMax(2)], ["uz"], -1);
            
            obj.x=ones(1,obj.analysis.getTotalElemsNumber());
            obj.result_number=13;
        end

        function setTempVars(obj,x)            
            for k=1:obj.nTempVars
                ni = find((round((obj.mesh.nodes(:,3)-obj.zCoords(k))*10000)/10000)==0);
                obj.fe.props.ndT(ni)=x(k);
            end
            obj.analysis.loadElementsThermal(obj.allganElemsSelector,obj.alphaT*obj.dT);
            %obj.analysis.loadElementsThermal(obj.allganTopElemsSelector,obj.alphaT*0.4*obj.dT);

        end

        function generateMesh( obj, ganTh, alGanTh, notchWidth, relNotchDepth, relRoutndNotchDepth )
            obj.ganTh=ganTh;
            obj.alGanTh=alGanTh;
            obj.notchWidth=notchWidth;
          
            
            % interface thickness
            obj.intTh = 0.025;
            
            %depth of the notch
            notchDepth=relNotchDepth*alGanTh; 
            
            %depth of the part ow the notch which is rounded
            roundNotchDepth=relRoutndNotchDepth*notchDepth; 
            
            %depth of the part ow the notch which is straight  
            straightNotchDepth=notchDepth-roundNotchDepth;

            % middle point of rounded noth relative height
            middleRoundedHeight = 0.15*roundNotchDepth;

            % middle point distance from noth bottom (middle of the notch)
            middleRoundedPlaneDistance = 0.3*notchWidth;

             % x - width of the tile
            obj.tileXWidth=30;
        
            % y - width of the tile
            obj.tileYWidth=20;


            % coordinates of shape points

            xr1=middleRoundedPlaneDistance;
            x1=notchWidth/2;
            x2=x1+obj.tileXWidth;
            x3=obj.tileXWidth+obj.notchWidth;
            xr2=x3-middleRoundedPlaneDistance;

            yr1=middleRoundedPlaneDistance;
            y1=notchWidth/2;
            y2=y1+obj.tileYWidth;
            y3=obj.tileYWidth+obj.notchWidth;
            yr2=y3-middleRoundedPlaneDistance;

            obj.xt=[ x3/2 y3/2 0 ];
            obj.xn=[ x3 y3 0 ];


            % top of gan layer
            zGt = ganTh;
            
            % top of interface layer
            zIt = zGt; +obj.intTh;
            
            % bottom of noth
            zNb = zIt+obj.alGanTh-notchDepth;
            
            %rounded middle point
            zRm = zNb+middleRoundedHeight;

            %rounded notch top point
            zRt = zNb+roundNotchDepth;

            %top of the structure
            zTop = obj.ganTh+obj.intTh+obj.alGanTh;
            

            % Mesh resolutions
            % number of tiles in the x direction
            xtiles=4;
            
            % number of tiles in the y direction
            ytiles=4; 
            
            % FE x - division of the tile
            ncx=8;
            
            % FE y - division of the tile
            ncy=8;
            
            % depth FE division of the GaN layer
            ngan=2;
            
            % FE depth division of the rouned part of the notch
            nround=1;
            
            % FE depth division of the straight part of the notch
            nstr=1;
            
            % FE  width division of the half notch (second half is symetrical if inside)
            nnotch=2;
            
            % value of the chemistry imposed
            chemistry=0.1;
            
          
            % One tile model generation
            ShapeFn8 = ShapeFunctionL8;
            ShapeFn27 = ShapeFunctionL27;

            mesh = Mesh();
            mesh2 = Mesh();
            
            %gan layer objects
            
            ganTileGeom = [  x1 y1 0;    x2 y1   0;  x1 y2   0;  x2 y2   0; ...
                             x1 y1 zGt;  x2 y1 zGt;  x1 y2 zGt;  x2 y2 zGt];

            ganNotch1 = [ 0 0 0;   x3 0 0;   x1 y1 0;   x2 y1 0; ...
                          0 0 zGt; x3 0 zGt; x1 y1 zGt; x2 y1 zGt ];

            ganNotch2 = [ 0 y3 0;   0 0 0;   x1 y2 0;   x1 y1 0; ...
                          0 y3 zGt; 0 0 zGt; x1 y2 zGt; x1 y1 zGt];
             

            thTileGeom = [  x1 y1 zGt;    x2 y1   zGt;  x1 y2   zGt;  x2 y2   zGt; ...
                            x1 y1 zIt;  x2 y1 zIt;  x1 y2 zIt;  x2 y2 zIt];

            thNotch1 = [ 0 0 zGt;   x3 0 zGt;   x1 y1 zGt;   x2 y1 zGt; ...
                          0 0 zIt; x3 0 zIt; x1 y1 zIt; x2 y1 zIt ];

            thNotch2 = [ 0 y3 zGt;   0 0 zGt;   x1 y2 zGt;   x1 y1 zGt; ...
                          0 y3 zIt; 0 0 zIt; x1 y2 zIt; x1 y1 zIt];

            allGanRound = [ x1 y1 zIt;   x2 y1 zIt;   x1 y2 zIt;   x2 y2 zIt; ...
                            x1 y1 zRt;   x2 y1 zRt;   x1 y2 zRt;   x2 y2 zRt];
            
            allGanStright = [ x1 y1 zRt;   x2 y1 zRt;   x1 y2 zRt;   x2 y2 zRt; ...
                            x1 y1 zTop;   x2 y1 zTop;   x1 y2 zTop;   x2 y2 zTop];
            
            zrm1=(zNb+zIt)/2;
            zrm2=(zRm+zIt)/2;
            zrm3=(zRt+zIt)/2;

            rNotch1 = [ 0   0   zIt;  x3/2 0   zIt; x3  0   zIt; ...
                        x1/2 y1/2 zIt;  x3/2 y1/2 zIt; x3-x1/2 y1/2 zIt; ...
                        x1  y1  zIt;  x3/2 y1  zIt; x2 y1   zIt; ...
                        0   0   zrm1;  x3/2 0   zrm1; x3  0   zrm1; ...
                        xr1 yr1 zrm2;  x3/2 yr1 zrm2; xr2 yr1 zrm2; ...
                        x1  y1  zrm3;  x3/2 y1  zrm3; x2 y1   zrm3; ...
                        0   0   zNb;  x3/2 0   zNb; x3  0   zNb; ...
                        xr1 yr1 zRm;  x3/2 yr1 zRm; xr2 yr1 zRm; ...
                        x1  y1  zRt;  x3/2 y1  zRt; x2 y1   zRt];

            rNotch2 = [ 0    y3      zIt; 0    y3/2 zIt; 0    0    zIt; ...
                        x1/2 y3-y1/2 zIt; x1/2 y3/2 zIt; x1/2 y1/2 zIt; ...
                        x1   y2      zIt; x1   y3/2 zIt; x1 y1 zIt; ...
                        0    y3      zrm1; 0    y3/2 zrm1; 0    0    zrm1; ...
                        xr1  yr2     zrm2; xr1  y3/2 zrm2; xr1 yr1   zrm2; ...
                        x1   y2      zrm3; x1   y3/2 zrm3; x1 y1     zrm3; ...
                        0    y3      zNb; 0     y3/2 zNb; 0    0     zNb; ...
                        xr1  yr2     zRm; xr1   y3/2 zRm; xr1 yr1    zRm; ...
                        x1   y2      zRt; x1    y3/2 zRt; x1  y1     zRt; ...

                ];
            
            
            mesh.addShapedMesh3D( ShapeFn8, ganNotch1, [ncx,nnotch,ngan], ShapeFn27.localNodes );
            mesh.addShapedMesh3D( ShapeFn8, ganNotch2, [ncy,nnotch,ngan], ShapeFn27.localNodes );
            mesh.addShapedMesh3D( ShapeFn8, thNotch1, [ncx,nnotch,1], ShapeFn27.localNodes );
            mesh.addShapedMesh3D( ShapeFn8, thNotch2, [ncy,nnotch,1], ShapeFn27.localNodes );
            mesh.addShapedMesh3D( ShapeFn27,rNotch1,  [ncx,nnotch,nround], ShapeFn27.localNodes );
            mesh.addShapedMesh3D( ShapeFn27,rNotch2,  [ncy,nnotch,nround], ShapeFn27.localNodes );
            mesh.duplicateTransformedMeshDeg3D( [x3/2  y3/2 ], 180, [0 0 0] );
            mesh.addShapedMesh3D( ShapeFn8, ganTileGeom, [ncx,ncy,ngan], ShapeFn27.localNodes );
            %mesh.addShapedMesh3D( ShapeFn8, thTileGeom,  [ncx,ncy,1], ShapeFn27.localNodes );
            mesh.addShapedMesh3D( ShapeFn8, allGanRound, [ncx,ncy,nround], ShapeFn27.localNodes );
            mesh.addShapedMesh3D( ShapeFn8, allGanStright, [ncx,ncy,nstr], ShapeFn27.localNodes );
            
            % copy tiles in x direction
            mesh.array(1,xtiles-1);
            
            % copy tiles in y direction
            mesh.array(2,ytiles-1);
            obj.mesh=mesh;
        end

        function s = getNodalStress(obj,result_node,result_number)
            sHM=obj.fe.results.nodal(result_node,result_number);
        end

        function computeZNodalCoords(obj)
            zCornersPoints=[    obj.mesh.nodes(obj.mesh.elems(:,1),:); ...
                                obj.mesh.nodes(obj.mesh.elems(:,3),:); ...
                                obj.mesh.nodes(obj.mesh.elems(:,7),:); ...
                                obj.mesh.nodes(obj.mesh.elems(:,9),:); ...
                                obj.mesh.nodes(obj.mesh.elems(:,19),:); ...
                                obj.mesh.nodes(obj.mesh.elems(:,21),:); ...
                                obj.mesh.nodes(obj.mesh.elems(:,25),:); ...
                                obj.mesh.nodes(obj.mesh.elems(:,27),:); ];
              
            obj.zCornersCoords = sort(unique(round(zCornersPoints(:,3)*obj.zTol)/obj.zTol),'descend');
        end

        function plotZCoordsPoints(obj)
            n=size(obj.zCornersCoords,1);
            x=zeros(n,1);
            p=plot3(x,x,obj.zCornersCoords,'.');
            p.Color = "red";
        end

        function [stressObj, sx1, sy1, sx2, sy2]=computeStressObjective(obj)
            n1 = obj.mesh.findClosestNode(obj.xt);
            n2 = obj.mesh.findClosestNode(obj.xn);
            sx1 = obj.fe.results.nodal.all(n1,7);
            sy1 = obj.fe.results.nodal.all(n1,8);
            sx2 = obj.fe.results.nodal.all(n2,7);
            sy2 = obj.fe.results.nodal.all(n2,8);
            stressObj= -( sx1+sy1 - (sx2+sy2) );
            %stressObj=(sx2+sy2)-(sx1+sy1);
        end

        % function o = computeObjectiveValue(obj,x)
        %     ganTh=h0*(1-x(1));
        %     alGanTh = h0*x(1);
        %     w=ganTh*x(2);
        %     r=h0*x(3);
        %     generateMesh( obj, ganTh, alGanTh, x(3), x(2) )
        % end

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
        
            fprintf(myfile,"\n EDIS\n");
            fprintf(myfile,"  gap 0.0001\n")
            for k=1:size(obj.zCornersCoords,1)
                if k<=obj.nTempVars
                    fprintf(myfile,"  3  %5.3f  0  0  0  %1.2f\n", obj.zCoords(k), chemistry(k) ); 
                else
                    fprintf(myfile,"  3  %5.3f  0  0  0  %1.2f\n", obj.zCoords(k), 0.0 ); 
                end
            end
                  
          meshMax=max(mesh.nodes);
          n1 = mesh.findClosestNode(  [0 0 0] );
          n2 = mesh.findClosestNode(  [meshMax(1) 0 0] );
          n3 = mesh.findClosestNode(  [meshMax(1) meshMax(2) ganTh] );
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

