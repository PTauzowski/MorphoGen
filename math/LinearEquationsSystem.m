classdef LinearEquationsSystem < handle
        
    properties
        I,J,supports,freedofs,supdofs,newdofs;
        iK11,iK12,iK21,iK22;
    end
          
    methods
        
        function obj = LinearEquationsSystem(I,J,supports)
           obj.I=I;
           obj.J=J;
           obj.supports = supports;
           obj.freedofs = find(not(obj.supports));
           obj.supdofs = find(obj.supports);
           ndofs = size(obj.freedofs,1);
           nsup = size(obj.supdofs,1);
           obj.newdofs = zeros(ndofs+nsup,1);
           obj.newdofs( obj.freedofs ) = (1:ndofs)';
           obj.newdofs( obj.supdofs ) = (1:nsup)';
           nI = obj.supports( obj.I );
           nJ = obj.supports( obj.J );
           obj.iK11 = not(nI) & not(nJ);
           obj.iK12 = not(nI) & nJ;
           obj.iK21 = nI & not(nJ);
           obj.iK22 = nI & nJ;
        end

        function q = solveClassical(obj,K,P)
             K11=sparse(obj.I,obj.J,K);
             maxK = max(max(K11));
             for k=1:max(size(obj.supdofs,1))
                K11( obj.supdofs(k), obj.supdofs(k) ) = maxK * 1.0e10;
             end
             q = K11 \ P;
        end
        function q = solve(obj,K,P)
             q=P;
             q(:)=0;
             dimfree  = size( obj.freedofs, 1 );
             q( obj.freedofs,:) =  sparse(obj.newdofs(obj.I(obj.iK11)),obj.newdofs(obj.J(obj.iK11)),K(obj.iK11),dimfree,dimfree) \ P( obj.freedofs,:);
        end
        function [qforms, lambdas] = solveEigenproblem(obj,K,Kg,num_eigenvalues)             
             dimfree  = size( obj.freedofs, 1 );
             [qforms, lambdas] = eigs( sparse(obj.newdofs(obj.I(obj.iK11)),obj.newdofs(obj.J(obj.iK11)),K( obj.iK11),dimfree,dimfree),...
                                       sparse(obj.newdofs(obj.I(obj.iK11)),obj.newdofs(obj.J(obj.iK11)),Kg(obj.iK11),dimfree,dimfree),...
                                       num_eigenvalues, 'smallestabs');
        end
        function [q, R, error] = solveR(obj,K,P)
            q=P;
            q(:)=0;
            R=P;
            R(:)=0;
            dimfree  = size( obj.freedofs, 1 );
            dimfixed = size( obj.supdofs, 1 );
            K11=sparse(obj.newdofs(obj.I(obj.iK11)),obj.newdofs(obj.J(obj.iK11)),K(obj.iK11),dimfree,dimfree);
            K21=sparse(obj.newdofs(obj.I(obj.iK21)),obj.newdofs(obj.J(obj.iK21)),K(obj.iK21),dimfixed,dimfree);
            q( obj.freedofs,:) = K11 \ P(obj.freedofs,:);
            R( obj.supdofs,:) = K21 * q( obj.freedofs,:);
            error = [ K11 * q(obj.freedofs) - P(obj.freedofs); K21 * q(obj.freedofs) - R ];
        end
        function q = solvePq(obj,K,P,q0)
             q=P;
             dimfree  = size( obj.freedofs, 1 );
             dimfixed = size( obj.supdofs, 1 );
             K11=sparse(obj.newdofs(obj.I(obj.iK11)),obj.newdofs(obj.J(obj.iK11)),K(obj.iK11),dimfree,dimfree);
             K12=sparse(obj.newdofs(obj.I(obj.iK12)),obj.newdofs(obj.J(obj.iK12)),K(obj.iK12),dimfree,dimfixed);
             q( obj.freedofs,:) = K11\(P( obj.freedofs,:)-K12*q0( obj.supdofs,:));
        end
        function [q, R, error] = solvePRq(obj,K,P,q0)
            q=P;
            q(:)=0;
            R=P;
            R(:)=0;
            dimfree  = size( obj.freedofs, 1 );
            dimfixed = size( obj.supdofs, 1 );
            obj.q0 = q0( obj.supdofs,:);
            obj.P = P( obj.freedofs,:); 
            K11=sparse(obj.newdofs(obj.I(obj.iK11)),obj.newdofs(obj.J(obj.iK11)),K(obj.iK11),dimfree,dimfree);
            K12=sparse(obj.newdofs(obj.I(obj.iK12)),obj.newdofs(obj.J(obj.iK12)),K(obj.iK12),dimfree,dimfixed);
            K21=sparse(obj.newdofs(obj.I(obj.iK21)),obj.newdofs(obj.J(obj.iK21)),K(obj.iK21),dimfixed,dimfree);
            K22=sparse(obj.newdofs(obj.I(obj.iK22)),obj.newdofs(obj.J(obj.iK22)),K(obj.iK22),dimfixed,dimfixed);
            q( obj.freedofs,:) = K11\(P( obj.freedofs )-K12*q0( obj.supdofs ));
            R( obj.supdofs,:) = K21 * q( obj.freedofs ) + K22 * q0( obj.supdofs );
            error = [ K11 * q( obj.freedofs, : ) + K12*q0( obj.supdofs, : ) - P( obj.freedofs, : ); ...
                      K21 * q( obj.freedofs, : ) + K22*q0( obj.supdofs, : ) - R( obj.supdofs, : ) ];
        end
       
    end
end

