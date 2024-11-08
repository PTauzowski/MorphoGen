classdef FrameOnElasticGroundModel
    
    properties
        Ecl, hcl, Ebm, hbm, b, lspan, hfloor, nspan, nfloor, Egr, nugr, lgr, hgr, resgr;
    end
    
    methods
        function obj = FrameOnElasticGroundModel( Ecl, hcl, Ebm, hbm, b, lspan, hfloor, nspan, nfloor, Egr, nugr, lgr, hgr, resgr, shapeFn2D )
            
                obj.Ecl=Ecl;
                obj.hcl=hcl;
                obj.Ebm=Ebm;
                obj.hbm=hbm;
                obj.b=b;
                obj.lspan=lspan;
                obj.hfloor=hfloor;
                obj.nspan=nspan;
                obj.nfloor=nfloor;
                obj.Egr=Egr;
                obj.nugr=nugr;
                obj.lgr=lgr;
                obj.hgr=hgr;
                obj.resgr=resgr;

                xp=[(lgr-nspan*lspan)/2 hgr];
       
                mesh=Mesh();
                frameFe = Frame2D( mesh.addHframe(nspan,lspan,nfloor,hfloor,xp) );
                frameFe.plot(mesh.nodes);
                elemsGround = mesh.addRectMeshArray2D(0,0,[xp(1) repelem(lspan,1,nspan) xp(1)], xp(2), resgr, shapeFn2D.pattern);

                elem2D = PlaneStressElem(shapeFn2D, elemsGround );
                elem2D.plot(mesh.nodes);
                
                nnodes  = size(mesh.nodes,1);
                nodesFr = false(nnodes,1);
                nodes2D = false(nnodes,1);
                nodesFr( frameFe.elems(:) ) = true;
                nodes2D( elem2D.elems(:) ) = true;
                find( nodes2D & nodesFr );

        end
        
        function outputArg = method1(obj,inputArg)
            
        end

    end
end

