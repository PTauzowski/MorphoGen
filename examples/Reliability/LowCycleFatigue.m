classdef LowCycleFatigue
    %LOWCYCLEFATIGUE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        E;
        sy;
        Cy;
        su;
        epD;
        m;
        sfi;
        S;
        s;  
        Dc;
        Nexp;
        Fmax;
        blocks;
        Fref;
    end
    
    methods
        function obj = LowCycleFatigue(fatigue_data)
            obj.E=fatigue_data.E;
            obj.sy=fatigue_data.sy;
            obj.Cy=fatigue_data.Cy;
            obj.su=fatigue_data.su;
            obj.epD=fatigue_data.epD;
            obj.m=fatigue_data.m;
            obj.sfi=fatigue_data.sfi;
            obj.S=fatigue_data.S;
            obj.s=fatigue_data.s;  
            obj.Dc=fatigue_data.Dc;
            obj.Nexp=fatigue_data.Nexp;
            obj.Fref = fatigue_data.Fref;
            %obj.Fmax = [35 40 45 50 55 60 65 70 75 80 85 90 95] / 1000;
            obj.Fmax = [50 50 50 50 50 50 50 50 50 50 50 50 50] / 1000;
            obj.blocks = [ 38000 100 100 100 100 100 100 100 100 100 100 100 ];
        end

        function ncycles = nCyclesMulti(obj,seqref)
            ncycles = zeros( 1, max(size(seqref)));
            for k=1:max(size(seqref))
                ncycles(k)=obj.nCycles(seqref(k));
            end
        end
        
        function ncycles = nCycles(obj,seqref)
            alphamax = obj.Fmax / obj.Fref;
            Eseelas = (alphamax*seqref).^2;
            smax=sqrt( (Eseelas+obj.E/obj.Cy*obj.sy^2)/(1+obj.E/obj.Cy));
            ds=sqrt((Eseelas+4*obj.E/obj.Cy*obj.sy^2)/(1+obj.E/obj.Cy));
            inds=find(ds>2*obj.sy);
            last = size(inds,2);
            if last==0
                 ncycles=1.0E6;
                return;
            end
            if size(inds,2)>1
                tblocks = obj.blocks(inds(1:end-1));
                sblocks = tblocks;
                sm = sblocks(1);
                for k=2:size(sblocks,2)
                    sm=sm+sblocks(k);
                    sblocks(k)=sm;
                end
                dep = (ds(inds)-2*obj.sy)/obj.Cy;
                dpdN = 2*(ds(inds)-2*obj.sy)/obj.Cy;
                %p=sum( obj.blocks(inds(1):last-1).*dpdN( inds(1):last-1) ) ;
                pD = obj.epD*( (obj.su - obj.sfi)./(ds/2-obj.sfi) ).^obj.m;
                Nd = pD(inds)./dep/2;
                sm=Nd(1);
                ii=0;
                for k=2:size(sblocks,2)
                    if ( sm<sblocks(k) )
                        ii = k-1;
                        Ndb = sblocks(k)-sm;
                        break;
                    end
                    sm=sm+Nd(k);
                end

                dDdN = (smax(inds).^(2*obj.s)+(ds(inds)-smax(inds)).^(2*obj.s))/(2*(2*obj.E*obj.S)^obj.s).*dpdN;
                D=0;
                if ii>0 
                    nblocks = obj.blocks(ii(1):size(obj.blocks,2));
                    nblocks(1)=Ndb;
                    D=sum( nblocks(inds(ii):last-1).*dDdN(inds(ii):last-1) );
                end
                Nr13=(obj.Dc-D)/dDdN(last);
                %pr = p+Nr13*dpdN(last);
                ncycles=sum(obj.blocks(inds(1:last-1)))+Nr13;
                42890-obj.Nexp;
                %abs(obj.fd.Nexp-ncycles)
            else
                dep = (ds(inds)-2*obj.sy)/obj.Cy;
                dpdN = 2*(ds(inds)-2*obj.sy)/obj.Cy;
                pD = obj.epD*( (obj.su - obj.sfi)./(ds/2-obj.sfi) ).^obj.m;
                dDdN = (smax(inds).^(2*obj.s)+(ds(inds)-smax(inds)).^(2*obj.s))/(2*(2*obj.E*obj.S)^obj.s).*dpdN;
                Nd = pD(inds)./dep/2;
                Nr =(obj.Dc)./dDdN(1);
                ncycles=Nd+Nr;
            end
        end
    end
end

