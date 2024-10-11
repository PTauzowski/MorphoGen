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
                baseRes = min( [hgr, lspan, xp(1)] );

                mesh=Mesh();
                frameElems = mesh.addHframe(nspan,lspan,nfloor,hfloor,xp);
                frameFe = Frame2D(frameElems);
                frameFe.plot(mesh.nodes);
                meshgr=Mesh();
                elemsgr = meshgr.addRectMesh2D(0,0,xp(1),hgr,round(xp(1)/baseRes*resgr),round(hgr/baseRes*resgr),shapeFn2D.pattern);
                xs=xp(1);
                for k=1:nspan
                    elemsgr = [elemsgr; meshgr.addRectMesh2D(xs,0,lspan,hgr,round(lspan/baseRes*resgr),round(hgr/baseRes*resgr),shapeFn2D.pattern)];
                    xs=xs+lspan;
                end
                elemsgr = [elemsgr; meshgr.addRectMesh2D(xs,0,xp(1),hgr,round(xp(1)/baseRes*resgr),round(hgr/baseRes*resgr),shapeFn2D.pattern)];
                elemsgr= mesh.merge(meshgr.nodes,elemsgr);
                
                elem2D = PlaneStressElem(shapeFn2D,elemsgr);
                elem2D.plot(mesh.nodes);

                


        end
        
        function outputArg = method1(obj,inputArg)
            
        end
    end
end

