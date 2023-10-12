classdef LinearEquationsSystemTr2D < handle
        
    properties
        I,J,supports,rotations,freedofs,supdofs,newdofs;
        iK11,iK12,iK21,iK22;
    end
          
    methods
        
        function obj = LinearEquationsSystemTr2D(I,J,supports,rotations)
           obj.I=I;
           obj.J=J;
           obj.supports = supports;
           obj.freedofs = find(not(obj.supports));
           obj.supdofs = find(obj.supports);
           obj.rotations = rotations;
           
        end

        function T = transformationMatrix(obj,K)
            dim = size(K,1);
            Td=sparse(1:dim,1:dim,ones(1,dim));
            T=sparse(1:dim,1:dim,ones(1,dim));
            irots = find(obj.rotations);
            for k=1:size(irots,1)
                alpha = obj.rotations(irots(k));
                ix=2*irots(k)-1;
                iy=2*irots(k);
                Ti=Td;
                Ti(ix,ix)=cos(alpha*pi/180);
                Ti(iy,iy)=cos(alpha*pi/180);
                Ti(ix,iy)=-sin(alpha*pi/180);
                Ti(iy,ix)=sin(alpha*pi/180);
                T=T*Ti;
            end
        end

        function q = solveClassical(obj,K,P)
             K11=sparse(obj.I,obj.J,K);
             T=obj.transformationMatrix(K11);
             K11r = T'*K11*T;
             K11=[];
             Pr = T'*P;
             maxK = max(max(K11r));
             for k=1:max(size(obj.supdofs,1))
                K11r( obj.supdofs(k), obj.supdofs(k) ) = maxK * 1.0e10;
             end
             qr = K11r \ Pr;
             q = T*qr;
        end
        function q = solve(obj,K,P)
             q=P;
             q(:)=0;
             Kt=sparse(obj.I,obj.J,K);
             T=obj.transformationMatrix(Kt);
             Ktr = T'*Kt*T;
             Kt=[];
             K11=Ktr(obj.freedofs,obj.freedofs);
             Ktr=[];
             P=T'*P;
             q(obj.freedofs,:) = K11 \ P( obj.freedofs, : );
             q = T*q;
        end
        function [q, R, error]  = solveR(obj,K,P)
            q=P;
            q(:)=0;
            R=P;
            R(:)=0;
            Kt=sparse(obj.I,obj.J,K);
            T=obj.transformationMatrix(Kt);
            Ktr = T'*Kt*T;
            Kt=[];

            Pr=T'*P;

            K11=Ktr(obj.freedofs,obj.freedofs);
            K21=Ktr(obj.supdofs,obj.freedofs);
            K22=Ktr(obj.supdofs,obj.supdofs);
            Ktr=[];

            q(obj.freedofs) = K11 \ Pr( obj.freedofs );
            R(obj.supdofs)  = K21 * q(obj.freedofs);
            error = [ K11 * q(obj.freedofs) - P(obj.freedofs); K21 * q(obj.freedofs) - R ];
            q=T*q;
            R=T*R;
        end
        function q = solvePq(obj,K,P,q0)
            q=P;
            q(:)=0;
            
            dim = dimfree+dimfixed;
            Kt=sparse(obj.I,obj.J,K);
            T=obj.transformationMatrix(Kt);
            Ktr = T'*Kt*T;
            Kt=[];

            Pr=T'*P;
            q0r=T'*q0;

            K11=Ktr(obj.freedofs,obj.freedofs);
            K12=Ktr(obj.freedofs,obj.supdofs);

            q(obj.freedofs,:) = K11 \ (Pr(obj.freedofs,:) - K12*q0r(obj.supdofs,:));
            q=T*q;
            
        end
        function [q, R, error] = solvePRq(obj,K,q0,P)
            q=P;
            q(:)=0;
            R=P;
            R(:)=0;
            dimfree  = size( obj.freedofs, 1 );
            dimfixed = size( obj.supdofs, 1 );
            dim = dimfree+dimfixed;
            
            Kt=sparse(obj.I,obj.J,K);
            T=obj.transformationMatrix(Kt);
            Ktr = T'*Kt*T;
            Kt=[];

            P=T'*P;
            q0=T'*q0;
            
            K11=Ktr(obj.freedofs,obj.freedofs);
            K12=Ktr(obj.freedofs,obj.supdofs);
            K21=Ktr(obj.supdofs,obj.freedofs);
            K22=Ktr(obj.supdofs,obj.supdofs);
            Ktr=[];

            q(obj.freedofs,:) = K11\(P(obj.freedofs,:)-K12*q0(obj.supdofs));
            R(obj.supdofs,:) = K21*q(obj.freedofs,:)+K22*q0(obj.supdofs);
            error = [ K11*q(obj.freedofs,:)+K12*q0(obj.supdofs)-P(obj.freedofd); K21*q(obj.freedofs,:)+K22*q0(obj.supdofs)-R(obj.supdofs) ];
            q=T*q;
            R=T*R;
        end
        
    end
end

