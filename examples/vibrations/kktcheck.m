function [residu,residunorm,residumax] = kktcheck(m,n,x,y,z,lam,xsi,eta,mu,zet,s,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d)
%KKTCHECK KKT residual check for MMA (standard).
ux1 = upp - x;  lx1 = x - low;
ux2 = ux1.^2;   lx2 = lx1.^2;

plam = p0 + P'*lam;
qlam = q0 + Q'*lam;
gvec = P*(1./ux1) + Q*(1./lx1);
f = a0*z + a'*y + (b + gvec);

residu1 = plam./ux2 - qlam./lx2 - xsi + eta;
residu2 = c + d.*y - mu - lam.*a;
residu3 = a0 - zet - a'*lam;
residu4 = f - y - z*ones(m,1);
residu5 = xsi.*(x-alfa);
residu6 = eta.*(beta-x);
residu7 = mu.*y;
residu8 = zet*z;
residu9 = lam.*s;

residu = [residu1;residu2;residu3;residu4;residu5;residu6;residu7;residu8;residu9];
residunorm = norm(residu,2);
residumax  = max(abs(residu));
end
