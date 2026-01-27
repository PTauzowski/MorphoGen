clear; clc; close all;
% Extracted from Huang et al. Appendix F/G (OCR) and cleaned.
% NOTE: If you encounter errors, compare with the PDF pages and fix remaining OCR artifacts.

% Helper: fast sparse (paper uses fsparse)
fsparse = @(i,j,s,sz) sparse(i,j,s,sz(1),sz(2));

% - ee
% PARAMETERS SETTING
DL = 8; DW = 4; DH = 4;
nelx = 80; nely = 40; nelz = 40;
vInt = [2.5 0.2 0.2 0 asin(1/sqrt(6)) asin(1/sqrt(5))];
volfrac = 0.2;
E = 2e11;nu = 0.3;rho = 7800;
mass = 1e6;
n_ord = 1;n_ord_all1l=6;
dgtO = 5;
scl_old = 1;
p = 6;lmd = 100;
iter = 1;maxiter = 500;
objVr5 = 1.0;
alpha = 1e-6; mid_num = 30; rank = 1; ini = 0.1; ep_old = 0.2;
loop = [];var_obj = [];Lamb_wi=[];Lamb_w2=[];
% hh ---- 7-7-7 o-oo %SEC 2): SETTING OF FE DISCRETIZATION
nEle = nelx*nely*nelz;
nNod = (nelx+1i)*(nely+1)*(nelz+1); nNodfc = (nelx+1)*(nely+1);
nDof = 3*nNod;
EL = DL/nelx; EW = DW/nely; EH = DH/nelz;
minSz = min([EL,EW,EH]);
[Ke] = Ke_tril(E,nu,EL,EW,EH);
KE(tril (ones (24))==1) = Ke'
KE = reshape (KE,24,24);KE = KE + KE' - diag(diag(KE));
[Me] = Me_tril(rho,EL,EW, EH);
ME(tril(ones (24))==1) = Me'
ME = reshape(ME,24,24);ME = ME + ME' - diag(diag(ME));
nodMat = int32(reshape(1:nNod,i+nelx,it+nely ,1+nelz));
edofVec = reshape (3*nodMat (1:nelx,1:nely ,1:nelz) ,nEle,1);
edofMat = edofVec + int32([3*nNodfc+[[3*nelx+[4 5 6 1 2 3]] -2 -1 0 1 2 3]...
% 3*nelx+[4 5 6 1 2 3] -2 -1 0 1 2 3]);
eleNodesID = edofMat(: ,[3:3:24]) ./3;
[sI,sII] = deal({]);
for j = 1:24
sI = cat(2,sI,j:24);
sII = cat(2,sII,repmat(j,1,24-j+1));
% end
[iK,jK] = deal(edofMat(:,sI)',edofMat(:,sII)');
TarO = sort([ik(:),jK(:)],2,'descend'); clear; ik jK
[x,y,z] = meshgrid(EL*[-nelx/2:nelx/2] ,EW*[-nely/2:nely/2] ,EH*[-nelz/2:nelz/2]) ;
LSgrid.x = permute(x,[2,1,3]); LSgrid.y = permute(y,[2,1,3]); LSgrid.z = permute(z
% 02) 1 31D);
volNod = sparse (double(eleNodesID(:)) ,1,1/8);
% #77 -c ttc cocoon %SEC 3): LOADS, DISPLACEMENT BOUNDARY CONDITIONS
io [jN,kN] = meshgrid(1:nely+1,1:nelz+1);
£ixNd = [1;(nelx+1)*(nely+1) ;nelx+1;nely*(nelx+1) +1];
fixDof = [3*fixNd(:)-2; 3*fixNd(:)-1; 3*fixNd(:)];
[jE,kE] = meshgrid(1:nely ,1:nelz);
fixEle = [1;nelx;nelx*nely;nelx*(nely-1) +1;
% (nelz-1) *nelx*nely+1;(nelz-1)*nelx*nelyt+nelx;(nelz-1)*nelx*nely+nelx*nely
% ;(nelz-1) *nelx*nely+nelx*(nely-1)+1];
freeDof = setdiff([1:nDof],fixDof);
massNd = (nely/2)*(nelx+1) + (nelx/2) + 1 + nelz*((nelx+1)*(nely+1));
massDof = 3*massNd-2:3*massNd;
massEle = (nelx/2) + (nely/2 - 1)*nelx + (nelz - 1)*nelx*nely;
% Ah morc corto coor ooo n%SEC 4): INITIAL SETTING OF COMPONENTS
x0=[kron(DL/4:DL/2:DL,ones(1,8)) kron(DL/4:DL/2:DL,ones(1,8))];
5s yO=[repmat (kron(DW/4:DW/2:DW,ones(1,4)),1,2) repmat (kron(DW/4:DW/2:DW, ones (1,4) )
% pil pAdI 8
zO=[DH/4*ones(1,16) 3*DH/4*ones(1,16)];
oo N = length(x0);

13 = repmat(vInt (3) ,1,N);

alp = repmat(vInt (4) ,1,N);

bet = repmat([1 -1 1 -1 1 -1 1 -1]*vInt(5),1,N/8);

gam = repmat([1 1 -1 -1 1 1 -1 -1]*vInt (6) ,1,N/8);

ov dd = [x0-DL/2; yO-DW/2; zO-DH/2; 11; 12; 13; alp; bet; gam];

nDsvb = length(dd(:));

nEhcp = nDsvb/N;

7o actComp = [1:N];

actDsvb = [1:nDsvb];

nNd = 0; PhiNd = [];

7s allPhi = [zeros(nNod,N) PhiNd];

% v4 fo ---- 0-2-2 %SEC 5): SETTING OF MMA

m= 13 c = 1000*ones(m,1); d = zeros(m,1);

aQ = 1; a = zeros(m,1);

xval = dd(:); xoldi = xval; xold2 = xval;

7s xmin = [-DL/2; -DW/2; -DH/2; minSz; minSz; minSz; -pi; -pi; -pil;

vo xmax = [DL/2; DW/2; DH/2; sqrt(DL*2+DW2+DH2)/2*[1; 1; 1]; pi; pi; pil;

xmin = repmat(xmin,N,1); xmax = repmat(xmax,N,1);

si low = xmin; upp = xmax;

% hh mr o onto %SEC 6): OPTIMIZATION LOOP

while objVr5 > 1e-4 && iter <= maxiter

Time_iteration=clock;

epsilon=initep_old -(1/ep_oldt+exp(-rank*(iter-mid_num)))~-1;

% omar tcc c rocco cco o-oo %LP 1): Generating TDFs and their derivatives

allPhiDrv = sparse (nNod,nDsvb) ;

for i = actComp

{allPhi,allPhiDrv,xval,actComp,actDsvb] = ...

% rN calc_Phi(allPhi ,allPhiDrv ,xval,i,LSgrid,p,nEhcp,epsilon,actComp,actDsvb
% »minSz);

% end

allPhiAct = [allPhi(:,actComp) PhiNd];

temp = exp(lmd*allPhiAct);

Phimax = max(-1e3,log(sum(temp ,2))/lmd) ;

allPhiDrvAct = allPhiDrv(:,actDsvb);

Phimaxdphi = kron(temp(:,1:length(actComp))./(sum(temp ,2)+eps) ,ones(1,nEhcp)) ;

PhimaxDrvAct = Phimaxdphi.*allPhiDrvAct;

% ocr c tcc cc cc ccc cco-o%LP 2): Plotting current design

% figure (1);clf

h = patch(isosurface(x,y,z,permute(reshape(Phimax, nelxt+i,nely+1,nelz+1)
»[2,1,3]) ,0));

hi = patch(isocaps(x,y,z,permute(reshape(Phimax, nelxt+i,nely+i1,nelz+1) ,[2,1,3])
20));

% set(h,'FaceColor','red','EdgeColor','none','facealpha',1); set(hi,'FaceColor','
% interp','EdgeColor','none');

% colormap([1 0 0]);

% isonormals(x,y,z,permute(reshape(Phimax, nelxti,nelyt+1,nelz+1) ,[2,1,3]),h);
% lighting flat;

% view(3) ,xlabel('x'),ylabel('y'),zlabel('z'), axis image;

% axis ([-DL/2,DL/2,-DW/2,DW/2,-DH/2,DH/2]) ; light ; pause (1e-1)

% drawnow

% omar tcc crc coco o o-oo %LP 3): Finite element analysis

H = Heaviside (Phimax ,alpha,epsilon) ;

den = sum(H(eleNodesID) ,2)/8;

Psil = zeros(nDof ,1);Psi2 = zeros(nDof ,1);Psi3 = zeros(nDof ,1);

bw_den=reshape (den ,nelx ,nely ,nelz) >alpha; bw_den_sol=reshape (bwlabeln(bw_den),
% nelx*nely*nelz,1);

bw_fix=unique (bw_den_sol (fixEle)) ;bw_mass=nonzeros (unique (bw_den_sol(massEle) ))
;

% clear; bw_den

% if sum(ismember (bw_fix,bw_mass))>0

struct=1;

eleLft=find(bw_den_sol==bw_mass) ;

denSld=zeros(size(den));

denSld(eleLft)=den(eleLft) ;

12( edofMatLft = edofMat(eleLft,:);

freedofLft = setdiff(edofMatLft ,fixDof);

Tar = sort ([iki(:),jki(:)],2,'descend'); clear; iki jK1

sK = reshape (Ke(:)*denSld(eleLft)',length(Ke)*length(eleLft) ,1);

K = fsparse(Iar(:,1),Iar(:,2),sk,[nDof ,nDof]); K = K + K' - diag(diag(K
% Ys

K = K + fsparse([1:nDof] ,[1:nDof],eps*ones(1,nDof) ,[nDof ,nDof]);

sM = reshape (Me(:)*denSld(eleLft)',length(Me)*length(eleLft) ,1);

M = fsparse(Iar(:,1),Iar(:,2),sM,[nDof ,nDof]); M=M + M' - diag(diag(M

M(massDof ,massDof)=M(massDof ,massDof)+mass*eye (length (massDof )) ;

M = M + fsparse([1:nDof],[1:nDof],eps*ones(1,nDof) ,[nDof ,nDof]);

Time_FEA=clock;

[PsiAll ,LambAll]=eigs(K(freedofLft ,freedofLft) ,M(freedofLft ,freedofLft),
% n_ord_all,'sm');

Lambi=LambAll1(n_ord ,n_ord);

Psil(freedofLft)=PsiAll(:,n_ord);

Psii=Psii/sqrt(Psil.'*M*Psil) ;

% else

% disp('WARNING!!! NO loading path is founded!!!');

sK=reshape (Ke(:)*den(:)', length (Ke) *nEle ,1);

K=fsparse(Iar0(:,1),Iar0(:,2),sK,[nDof ,nDof]) ;K=K+K'-diag(diag(K));

sM=reshape (Me(:)*den(:)',length(Me) *nEle ,1);

M=fsparse(TIar0(:,1),Iar0(:,2),sM,[nDof ,nDof]) ;M=M+M'-diag (diag (M));

M(massDof ,massDof)=M(massDof ,massDof )+mass*eye (length (massDof )) ;

Time_FEA=clock;

[PsiAll ,LambAll]=eigs(K(freeDof ,freeDof) ,M(freeDof ,freeDof) ,n_ord_all,'sm')

Lambi=LambAll (n_ord ,n_ord) ;

Psii(freeDof)=PsiAll(:,n_ord);

Psi1=Psil/sqrt(Psil.'*M*Psil);

% end

OBJ (iter)=Lambi;

fval=sum(den) *EW*EH*EL/(DW*DH*DL)-volfrac; CONS(iter) = fval + volfrac;

% eo matt ttt occ cca %LP 4): Sensitivity analysis

df0dx=zeros(1,nDsvb) ;dfdx=zeros(1,nDsvb) ;

delta_H=3*(1-alpha) /(4*epsilon) *(1-Phimax.2/(epsilon"2));

delta_H (abs (Phimax) >epsilon) =0;

if struct == 1

Cai_wi,a2_wiJ=calc_fOval(K,M,freedofLft ,Psii ,nDof ,Lambi);

% else

{ait_wi,a2_wi]J=calc_fOval(K,M,freeDof ,Psi1,nDof ,Lamb1);

% end

energy=sum((at_wi(edofMat) *(KE-Lamb1.*ME)) .*Psil(edofMat) -0.5*(a2_w1*Psii1(
% edofMat)*ME) .*Psil(edofMat) ,2) ;

sEner=energy*ones (1,8) /8;

engyNod=sparse (double (eleNodesID(:)),1,sEner(:));

df0dx (actDsvb)=(engyNod.*delta_H) '*PhimaxDrvAct;

dfdx(actDsvb)=(volNod.*delta_H) '*PhimaxDrvAct *EW*EH*EL/(DW*DH*DL) ;

dgt=dgt0-floor(1log10([max (abs (df0dx(:))) max(abs(dfdx(:)))]));

% if iter>1 && OBJ(iter)-OBJ(iter-1)>0

scl=max (abs (df0dx));scl_old=scl;

% else

scl=scl_old;

% end

f0val=Lamb1/scl;

df0dx=round (df0dx*10 dgt (1) )/10dgt (1)/scl;

dfdx=round (dfdx*10dgt (2) )/10 dgt (2);

(xmma,~,~,~,~,~,~,7,7,low,upp] = mmasub(m,nDsvb,iter,xval(:),...

% xmin ,xmax ,xoldi1,xold2,fOval ,df0dx ,fval ,dfdx,low,upp,a0,a,c,d);

xold2=xoldi;xoldi=xval;xval=xmma;

if iter>=5 && fval/volfrac<le-4

objVr5=abs (max (abs (OBJ (iter-4:iter)-mean(OBJ(iter-4:iter))))/mean(OBJ(iter
% -4:iter)));

% end

% disp(['? It.:  sprintf('?%4i\t',iter) ' Obj.: ' sprintf('%6.3f\t',f0val*scl) '
% Vol.: " ...

% sprintf('%6.4f\t',fval) ' Ch.: ' sprintf('%6.4f\t',objVr5) 'Time_iteration
% :?,num2str(etime(clock,Time_iteration))]);

% figure (2);

var_obj=[var_obj Lamb1];loop=[loop iter];

% yyaxis left

% plot (loop,var_obj ,'DisplayName','wi')

% xlabel('Step')

% ylabel(' Frequency ')

% yyaxis right

% plot (loop ,CONS ,'DisplayName','vol')

% ylabel('Vol')

% ylim(f[o 1])

% legend('show','location','southeast ')

% drawnow

iter=itert+1;

% end

os function [al,a2]=calc_fOval(K,M,Dof ,Psi,nDof ,Lamb)

L11=K(Dof ,Dof)-Lamb*M(Dof ,Dof) ;

L12=-M(Dof ,Dof)*Psi(Dof);L21=L12.';

LL = [L11,L12;

% L21,0];

vec_a=LL\[zeros(size(Dof ,1) ,1);1];

al=zeros(nDof ,1) ;a1(Dof)=vec_a(1:end-1) ;a2=vec_a(end);

% end

function [allPhi,allPhidrv,xval,actComp,actDsvb] = ...

% calc_Phi(allPhi ,allPhidrv ,xval,i,LSgrid,p,nEhcp,epsilon,actComp,actDsvb,minSz)

di = xval((i-1)*nEhcp+1:i*nEhcp) ;

xO = di(1); yO = di(2); zO = di(3); 11 = di(4) + eps; 12 = di(5) + eps; 13 = di
% (6) + eps;

sa = sin(di(7)); sb = sin(di(8)); sg = sin(di(9));

ca = cos(di(7)); cb = cos(di(8)); cg = cos(di(9));

R = [cb*cg cb*sg -sb; sa*sb*cg-ca*sg sa*sb*sgtca*cg sa¥cb;...

% ca*sb*cg+sa*sg ca*sb*sg-sa*cg ca*cb];

xyzLc = [LSgrid.x(:)-x0 LSgrid.y(:)-yO LSgrid.z(:)-z0]; % local
% coordinates of all nodes

xyz = xyzLc*R';

x1 = xyz(:,1) + eps; yl = xyz(:,2) + eps; zi = xyz(:,3) + eps;

temp = (x1/11).*p + (y1/12).*p + (21/13)."p;

allPhi(:,i) = 1 - temp.*(1/p); % TDF of i-th
% component

% if di(5)/minSz < 1.01

% disp(['The ' sprintf('%i',i) '-th component is too small! DELETE it!!!?]);

allPhi(:,i) = -1e3;

xval((i-1)*nEhcp+[4:6]) = 0;

actComp = setdiff(actComp,i);

actDsvb = setdiff (actDsvb ,nEhcp*i-nEhcpt1: nEhcp*i) ;

% elseif di(6)/minSz < 1.01 && di(4)/minSz < 1.01

% disp(['The ' sprintf('%i',i) '-th component is too small! DELETE it!!!']);

allPhi(:,i) = -1e3;

xval((i-1)*nEhcp+[4:6]) = 0;

actComp = setdiff(actComp,i);

actDsvb = setdiff(actDsvb ,nEhcp*i-nEhcp+1:nEhcp*i) ; elseif min(abs(
allPhi(:,i))) >= epsilon

% disp(['?The ' sprintf('%i',i) '-th component is too small! DELETE it!!!']);

allPhi(:,i) = -1e3;

xval((i-1)*nEhcp+[4:6]) = 0;

actComp = setdiff(actComp,i);

actDsvb = setdiff(actDsvb ,nEhcp*i-nEhcp+1:nEhcp*i) ;

% else

Ra = [0 0 0; R(3,1) R(3,2) R(3,3); -R(2,1) -R(2,2) -R(2,3)];

Rb = [-sb*cg -sb*sg -cb; sa*xcb*cg sa*cb*sg -sa*sb; ca*cb*cg ca*cb*sg -ca*sb
1;

Rg = [-cb*sg cb*cg 0; -sa¥*sb*sg-ca*cg sa*sb*cg-ca*sg 0; -ca*sb*sgtsa*cg ca*
% sb*cg+sa*sg 0];

dxi = [[-R(1,:) 0.0 0.0 0.0]+0.0*x1 xyzLc*[Ra(1,:); Rb(1,:); Rg(1,:)]7]; %
% variation of x'

% variation of y'

dzi = [[-R(3,:) 0.0 0.0 0.0]+0.0*y1 xyzLc*[Ra(3,:); Rb(3,:); Rg(3,:)]]; %
% variation of z'

temp1 = -temp.*(1/p-1).*(x1/11).*(p-1);

temp2 = -temp.*(1/p-1).*(y1/12) .*(p-1);

temp3 = -temp.~(1/p-1) .*(z1/13).*(p-1);

dpdxi = tempi/li; dpdy1l = temp2/12; dpdzi = temp3/13;

dpdli = -tempi.*(x1/112); dpdl2 = -temp2.*(y1/122); dpd1l3 = -temp3.¥*(z1/
1372);

dpdL = 0.0*dx1; dpdL(:,4:6) = [dpdli dpd12 dpd13];

allPhidrv(:,nEhcp*(i-1)+1:nEhcp*i) = dpdL + dpdxi.*dx1 + dpdy1.*dy1 + dpdzi
% .*dz1;

% end

% end

function [H] = Heaviside(phi,alpha,epsilon)

H = 3*(1-alpha)/4*(phi/epsilon-phi.*3/(3*(epsilon)~3)) + (1+alpha)/2;

H(phi>epsilon) = 1;

H(phi<-epsilon) = alpha;

% end

function [Ke] = Ke_tril(E,nu,EL,EW, EH)

ti = EL*EW/36/EH; t2 = EL*EH/36/EW; t3 = EW*EH/36/EL;

pi = EL/48; p2 = EW/48; p3 = EH/48;

ki = nux[-4*(t1+t2+t3) 0 0 2*(-t1-t2+2*t3) 8*p3 8*p2 -t1+2*(t2+t3) O 4*p2 -2*(
% t1-2*t2+t3) -8*p3 0...

% 2*(2*t1-t2-t3) 0 -8*p2 2*(t1+t3)-t2 4*p3 O t1+t2+t3 0 O 2*(t1+t2)-t3 -4*p3
% -4*p2] +...

% [2*(t1+t2+2*t3) 2*p3 2*p2 tit+t2-4*t3 -2*p3 -2*p2 t1/2-t2-2*t3 -2*p3 -p2 ti
% -2*(t2-t3) 2*p3 p2...

% -2«(t1-t3)+t2 p3 2ep2 -tit+t2/2-2*t3 -p3 -2*p2 -(ti+t2)/2-t3 -p3 -p2 -(ti+
% t2)+t3 p3 p2];

k2 = nu*[-4*(t1i+t2+t3) 0 -8*p3 -2*(t1+t2-2*t3) 0 O -t1+2*(t2+t3) 4*p1 B*p3 -2*(
% t1-2*t2+t3) 8*p1 0...

% 2*(2*t1-t2-t3) -8*p1 -4*p3 2*(t1+t3)-t2 -4*p1 0 t1+t2+t3 0 4*p3 2*(t1+t2)-
% t3 O] +...

% [2*(t1+2*t2+t3) 2*p1 2p3 t1i+2*(t2-t3) pl -2*p3 t1/2-2*t2-t3 -p1 -2*p3 t1
% -4*t2+t3 -2*p1 p3 ...

% 27¢ -2*(t1-t2)+t3 2*p1 p3 -t1+t2-t3 pi -p3 -(t1+t3)/2-t2 -p1 -p3 -t1-2*t2+t3/2
% -2*p1];

k3 = nu*[-4*(t1i+t2+t3) -8*p2 0 -2*(t1+t2-2*t3) -4*p2 -4*p1 -t1+2*(t2+t3) 0 -8*
% pi -2*(t1-2*t2+t3) ...

% 8*p2 8*pi 2*(2Q*t1-t2-t3) O 4*pi 2*(t1+t3)-t2 O O t1i+t2+t3 4*p2 O 2*(t1+t2)
=S8]] & ooo

% [2*(2*t1+t2+t3) 2*p2 pi 2*(t1-t3)+t2 p2 pi ti-t2-t3 p2 2*pi 2*(t1-t2)+t3
% -2*p2 -2*p1 ...

% -4*t1t+t2+t3 -2*p2 -p1 -2*t1+t2/2-t3 -p2 -pi -t1-(t2+t3)/2 -p2 -2*p1 -2*t1-
% t2+t3/2);

Ke = E/((1+nu) *(1-2*nu))*([k1?;k2?;k3?;k1(1) ;k1(8) 5 k1(18) ;k1(10) ;k1(5) 5k1(21) 5 k1
% (7) ;k1(2) 5;k1(24); ...

% k1 (16) ;k1(23) ;k1(3) ;k1(13) ;k1(20) ;k1(6) 5 k1 (22) 5k1(17) 5k1(9) 5k1(19) ;k1(14) ;
% k1(12) ;k2(1);k2(2); ...

% k2(3) ;k2(10) ;k2 (11) ;k1(2) 5 k2(7) 5; k2(8) 5k1(17) 5k2(16) 5k2 (17) 5k1(20) ;k2(13) ;
% k2(14) ;k1(23) ;k2(22); ...

% k2 (23) ;k1(14) ;k2(19) ;k2(20) ;k3 (1) ;k1(21) 5k2 (14) ;k3 (10) ;k1(9) ;k2(17) ;k3(7);
% k1(3) ;k2(8);k3(16); ...

% ki (15) ;k2(11) ;k3 (13) ;k1(24) ;k2 (23) ;k3 (22) ;k1(12) ;k2(20) ;k3 (19) ;k1(1) ;k1(2)
% 3k1(18) ;k1(4);k1(5); ...

% 28( k1 (15) ;k1(19) ;k1(20) ;k1(12) ;k1(22) ;k1(23) ;k1(9) ;k1(13) ;k1(14) ;k1(6) 5 k1 (16)
% 3k1(17) ;k1(3);k2(1); ...

% k2 (23) ;k2(3) ;k2(4) ;k2(20) ;k1(20) ;k2(19) ;k2(5) ;k1(17) ;k2(22) ;k2(2) ;k1(14);
% k2(13) ;k2(11);k1(23); ...

% k2 (16) ;k2(8) ;k3 (1) ;k1(6) ;k2 (20) ;k3(4) ;k1(12) ;k2(5) ;k3 (19) ;k1(24) ;k2(2) 5; k3
% (22) ;k1(15) ;k2(14); ...

% k3 (13) ;k1(3) ;k2(17) ;k3 (16) 5 k1(1) 5 k1(8) 5 k1(3) 5 k1(22) 5k1 (17) 5k1(24) ;k1(19) ;
% k1(14) 5k1(21);k1(16); ...