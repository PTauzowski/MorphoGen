classdef (Abstract) FEAnalysis < handle
    
    properties
        felems, mesh, ndofs, supports, rotations, Pnodal, qnodal, Pfem, qfem;
        selTolerance;
    end
  
    methods
        function obj = FEAnalysis(felems, mesh)
            if iscell(felems)
                obj.felems = felems;
            else
                obj.felems = { felems };
            end
            obj.mesh = mesh;
            obj.ndofs = obj.felems{1}.ndofs;
            for k=max(size(obj.felems))
                obj.ndofs = union(obj.ndofs,obj.felems{k}.ndofs);
            end
            obj.selTolerance=1.0E-05;
            obj.Pnodal = zeros( size(mesh.nodes,1), size(obj.ndofs,2) );
            obj.Pfem = zeros( size(mesh.nodes,1) * size(obj.ndofs,2), 1 );
            obj.qnodal = zeros( size(mesh.nodes,1), size(obj.ndofs,2) );
            obj.supports = zeros( size(mesh.nodes,1), size(obj.ndofs,2) );
        end
        function dim = getTaskDim(obj)
            dim = size(obj.mesh.nodes,1)*size(obj.ndofs,2);
        end
        function ne = getTotalElemsNumber(obj)
            ne=0;
            for k=1:size(obj.felems,1)
                ne = ne + size( obj.felems{k}.elems, 1);
            end
        end
        function ei = getElemIndices(obj)
             ne=size(obj.felems,1);
             ei=cell(ne,1);
             offset=1;
             for k=1:ne  
                 ei{k}=offset:offset+size(obj.felems{k}.elems,1)-1;
                 offset=offset+size(obj.felems{k}.elems,1);
             end
        end
        function selems = selectElems(obj,selector)
             ne=size(obj.felems,1);
             tne=obj.getTotalElemsNumber();
             selems=[];
             for k=1:ne  
                 selems = find( ...
                            selector( ...
                                permute( ...
                                    mean( ...
                                        reshape( ...
                                            obj.mesh.nodes( ...
                                                obj.felems{1}.elems', ...
                                                : ...
                                                ), ...
                                            size(obj.felems{1}.elems,2), ...
                                            size(obj.felems{1}.elems,1), ...
                                            size(obj.mesh.nodes,2) ...
                                         ) ...
                                     ), ...
                                     [3,2,1])' ...
                                 ) ...
                               );
             end
        end
        function setRotations(obj, r)
           obj.rotations=r;
        end
        function setNodalSupportRotations(obj, selector, angle)
           if size(obj.rotations,1) == 0
               obj.rotations=zeros(size(obj.mesh.nodes,1),1);
           end
           inodes=find(selector.select(obj.mesh.nodes));
           %inodes=selector;
           obj.rotations(inodes)=angle;
        end
        function Fe = toFEMVector( obj, F )
                Fe = reshape(F',obj.getTaskDim(),1);
        end
        function F = fromFEMVector( obj, Fe )
            F=reshape(Fe,size(obj.ndofs,2),size(obj.mesh.nodes,1))';
        end
        function K = globalMatrixAggregation(obj, fname)
            K = [];
            for k=1:max(size(obj.felems))
                if ismethod(obj.felems{k},fname)
                    K = [ K obj.felems{k}.(fname)(obj.mesh.nodes) ];
                else
                    error("Class " + class(fe) + " or its predecessors not implements function :"+fname);
                end
            end
        end
        function K = globalSolutionDependendMatrixAggregation(obj, fname, q)
            K = [];
            for k=1:max(size(obj.felems))
                if ismethod(obj.felems{k},fname)
                    K = [ K obj.felems{k}.(fname)(obj.mesh.nodes, q) ];
                else
                    error("Class " + class(fe) + " or its predecessors not implements function :"+fname);
                end
            end
        end
        function [I,J,V,Ksize] = globalMatrixIndices(obj)
            I=[];
            J=[];
            V=[];
            Ksize = 0;
            for k=1:max(size(obj.felems))
                [Ie,Je,Ve,Kesize] = obj.felems{k}.sparseMatrixAllocDataUniform( obj.ndofs );
                I = [ I reshape(Ie',[],1) ];
                J = [ J reshape(Je',[],1) ];
                V = [ V reshape(Ve',[],1) ];
                Ksize = Ksize + Kesize;
            end
        end
        function id = findDOFsIndices(obj,dofnames)
            [~,id,~] = intersect(obj.ndofs,dofnames);
        end
        function createNextRightHandSideVector(obj)
            obj.Pfem=[obj.Pfem obj.toFEMVector(obj.Pnodal) ];
            obj.Pnodal(:) = 0;
        end
        function setRightHandSideVectorAsCurrent(obj,n)
            obj.Pnodal = obj.fromFEMVector( obj.Pfem(:,n) );
        end
        function setCurrentLoadToRightHandSideVectors(obj,n)
            obj.Pfem(:,n) = obj.toFEMVector( obj.Pnodal );
        end
        function q = getSolution(obj,n)
            obj.qnodal = obj.fromFEMVector( obj.qfem(:,n) );
            q = obj.qnodal;
        end
        function P = getCurrentNodalLoad(obj)
            P=obj.Pnodal;
        end
        function P = setCurrentNodalLoad(obj,P)
            obj.Pnodal=P;
        end
        function prepareRHSVectors(obj)
            if any(obj.Pnodal(:)~=0) && size(obj.Pfem,2)>1
                obj.createNextRightHandSideVector();
            else
                setCurrentLoadToRightHandSideVectors(obj,size(obj.Pfem,2))
            end
        end
        function loadClosestNode(obj, x, dofnames, values )
           obj.Pnodal( obj.mesh.findClosestNode(x), obj.findDOFsIndices( dofnames ) ) = obj.Pnodal( obj.mesh.findClosestNode(x), obj.findDOFsIndices( dofnames ) ) + values;
        end
        function loadNodes(obj, nodesel, dofnames, values )
           snodes = find( nodesel.select( obj.nodes) );
           obj.Pnodal( snodes, obj.findDOFsIndices( dofnames ) ) =  obj.Pnodal( snodes, obj.findDOFsIndices( dofnames ) ) + repmat(values, size(snodes,1),1 );
        end
        function loadNodesByNumbers(obj, nodes, dofnames, values )
            snodes=nodes';
           obj.Pnodal( snodes, obj.findDOFsIndices( dofnames ) ) =  obj.Pnodal( snodes, obj.findDOFsIndices( dofnames ) ) + repmat(values, size(snodes,1),1 );
        end
        function elementLoadLineIntegral(obj, mode, edgesSel, dofnames, valueFn )
            di = obj.findDOFsIndices(dofnames);
            for k=1:max(size(obj.felems))
                fedges = obj.felems{k}.findEdges(obj.mesh.findNodes(edgesSel)); 
                obj.Pnodal = obj.Pnodal + obj.felems{k}.loadLineIntegral(mode, obj.mesh.nodes, fedges, dofnames, di, obj.Pnodal, valueFn);
            end
        end
        function elementLoadSurfaceIntegral(obj, mode, faceSel, dofnames, valueFn )
            di = obj.findDOFsIndices(dofnames);
            for k=1:max(size(obj.felems))
                ffaces = obj.felems{k}.findFaces(obj.mesh.findNodes(faceSel)); 
                obj.Pnodal = obj.Pnodal + obj.felems{k}.loadSurfaceIntegral(mode, obj.mesh.nodes, ffaces, dofnames, di, obj.Pnodal, valueFn);
            end
        end
        function fixClosestNode(obj, x, dofnames, values )
           node = obj.mesh.findClosestNode(x);
           obj.qnodal(node, obj.findDOFsIndices( dofnames ) ) = values;
           obj.supports(node, obj.findDOFsIndices( dofnames ) ) = ones( 1, size( dofnames, 2) );
        end
        function fixNodes(obj, selector, dofnames, varargin )
           inds = obj.findDOFsIndices( dofnames );
           snodes = find(selector.select(obj.mesh.nodes));
           if  nargin==4
               obj.q(snodes, inds) = repmat(varargin{1}, size(snodes,1),1 );
           end
           obj.supports(snodes, obj.findDOFsIndices( dofnames ) ) = ones( size(snodes,1), size(dofnames,2) );
        end
        function computeElementResults(obj,q,varargin)
            resnumber=0;
            nnodes=size(obj.mesh.nodes,1);
            for k=1:max(size(obj.felems))
                if ( nargin == 3 )
                    obj.felems{k}.computeResults( obj.mesh.nodes,obj.qnodal,varargin{1});
                    %obj.felems{k}.computeResults( obj.nodes,obj.fromFEMVector(q),varargin{1});
                else
                    obj.felems{k}.computeResults(obj.mesh.nodes,obj.qnodal);
                    %obj.felems{k}.computeResults(obj.nodes,obj.fromFEMVector(q));
                end
                GPresults=permute( obj.felems{k}.results.GPvalues,[3,1,2]);
                gpres = GPresults( 1,:,: );  
                resnumber=max(resnumber,size( gpres, 2 ));
            end
            nres = zeros( nnodes, resnumber );
            ires = zeros( nnodes, resnumber );
            for k=1:max(size(obj.felems))
                GPresults = permute( obj.felems{k}.results.GPvalues,[3,1,2]);
                gpres = GPresults( 1,:,: );  % gp x results
                el = obj.felems{k}.elems;
                sfv = obj.felems{k}.sf.getRecoveryMatrix();
                for k=1:size(el,1)
                  neres = sfv * GPresults(:,:,k);
                  nres( el( k, : ), : ) = nres( el( k, : ), : ) + neres;
                  ires( el( k, : ), : ) = ires( el( k, : ), : ) + 1;
                end
            end
            for k=1:max(size(obj.felems))
              obj.felems{k}.results.nodal = nres ./ ires;
            end
        end
        function printProblemInfo(obj)
            disp('');
            disp(' *** Matlab Object Oriented Finite Element Method *** ');
            disp(' *** P. Tauzowski (C) 2021 ************************* ');
            disp("Problem type                    :"+class(obj));
            disp("Number of nodes                 :"+num2str(size(obj.mesh.nodes,1)));
            disp("Number of degrees of freedom    :"+num2str(obj.getTaskDim()));
            disp("Number of element types         :"+num2str(max(size(obj.felems))));
            ne=0;
            for k=1:max(size(obj.felems))
                disp("     element "+class(obj.felems{k}));
                disp("     number of instances :"+num2str(size(obj.felems{k}.elems,1)));
                ne = ne + size(obj.felems{k}.elems,1);
            end
            disp("Total number of finite elements :"+num2str(ne));
            
        end
        function plotNodes(obj)
            hold on;
            daspect([1 1 1]);
            if  size( obj.mesh.nodes, 2) == 2 
                plot(obj.mesh.nodes(:,1), obj.mesh.nodes(:,2), 'o');
            elseif size( obj.mesh.nodes, 2) == 3
                 plot3(obj.mesh.nodes(:,1), obj.mesh.nodes(:,2), obj.mesh.nodes(:,3), 'o');
            end
                
        end
        function plotSelectedNodes(obj,sel)
            hold on;
            daspect([1 1 1]);
            sn = sel.select( obj.mesh.nodes );
            plot(obj.mesh.nodes(sn,1), obj.mesh.nodes(sn,2), 'o');
        end
        function plotCurrentLoad(obj)
              dim    = size(obj.mesh.nodes,2);
              dg     = norm( max(obj.mesh.nodes) - min(obj.mesh.nodes) );
              maxs = max( abs(min(min(obj.Pnodal))), abs(max(max(obj.Pnodal)) ) );
              xp = obj.mesh.nodes - obj.Pnodal ./ maxs * dg * 0.02;

              if dim == 2
                X = [xp(:,1) obj.mesh.nodes(:,1)]';
                Y = [xp(:,2) obj.mesh.nodes(:,2)]';
                plot( X, Y, 'm', 'LineWidth', 4 );
             elseif dim == 3
                X = [xp(:,1) obj.mesh.nodes(:,1)]';
                Y = [xp(:,2) obj.mesh.nodes(:,2)]';
                Z = [xp(:,3) obj.mesh.nodes(:,3)]';
                plot3( X, Y, Z, 'm', 'LineWidth', 4 );
             end
        end
        function clearCurrentLoad(obj)
                obj.Pnodal(:)=0;
        end
        function plotSupport(obj)
              dim    = size(obj.mesh.nodes,2);
              dg     = norm( max(obj.mesh.nodes) - min(obj.mesh.nodes) );
              xs =  obj.mesh.nodes;
              xps = xs;

              irots = find(obj.rotations);

              xp = xps - obj.supports * dg * 0.02;
              
              for k=1:size(irots)
                  alpha = obj.rotations(irots(k));
                  xp(irots(k),1) = xps(irots(k),1) - cos(alpha*pi/180) * dg * 0.02; 
                  xp(irots(k),2) = xps(irots(k),2) - sin(alpha*pi/180) * dg * 0.02; 
              end

             if dim == 2
                X = [xp(:,1) xs(:,1)]';
                Y = [xp(:,2) xs(:,2)]';
                plot( X, Y, 'b', 'LineWidth', 4 );
             elseif dim == 3
                X = [xp(:,1) xs(:,1)]';
                Y = [xp(:,2) xs(:,2)]';
                Z = [xp(:,3) xs(:,3)]';
                plot3( X, Y, Z, 'b', 'LineWidth', 4 );
             end
        end
        function plotMaps( obj, mapNames, scd )
            for mapName=mapNames(:)'
                figure, hold on, axis off;
                daspect([1 1 1]);
                dg     = norm( max(obj.mesh.nodes) - min(obj.mesh.nodes) );
                maxs = max( abs(min(min(obj.qnodal))), abs(max(max(obj.qnodal)) ) );
                for k=1:max(size(obj.felems))
                   obj.felems{k}.plotMap( obj.mesh.nodes, obj.qnodal, mapName, dg / maxs * scd );
                end
            end
        end
        function plotAllNodeNumbers(obj)
            text(obj.nodes(:,1),obj.nodes(:,2),num2str(1:size(obj.nodes,1)));
        end
        function plotSelectedNodeNumbers(obj,selector)
            ninds=find(selector.select(obj.nodes));
            snodes = obj.nodes(ninds,:);
            text(snodes(:,1),snodes(:,2),num2str(ninds));
        end
     end
        
end

