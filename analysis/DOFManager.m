classdef DOFManager < handle
   
    properties
        allDofTypes;
    end

     methods(Abstract)
         getIndices(obj,fElem);
     end
    
    methods

        function [I,J,V,Ksize] = globalMatrixIndices(obj,fElems)
            I=[];
            J=[];
            V=[];
            Ksize = 0;
            for k=1:max(size(fElems))
                [Ie,Je,Ve,Kesize] = obj.getIndices( fElems{k} );
                I = [ I reshape(Ie',[],1) ];
                J = [ J reshape(Je',[],1) ];
                V = [ V reshape(Ve',[],1) ];
                Ksize = Ksize + Kesize;
            end
        end


    end

end

