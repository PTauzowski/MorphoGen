clear; clc; close all;
% Extracted from Huang et al. Appendix F/G (OCR) and cleaned.
% NOTE: If you encounter errors, compare with the PDF pages and fix remaining OCR artifacts.

% Helper: fast sparse (paper uses fsparse)
fsparse = @(i,j,s,sz) sparse(i,j,s,sz(1),sz(2));

% PARAMETERS SETTING
DW=4;DH=1;h=0.01; % width/height/thickness of design domain
nelx=400;nely=100; % column/row number of finite elements
xInt=0.25; yInt=0.25; % distance between initial components' center in x/y
% direction
vInt=[0.4 0.04 0.04 pi/4]; % other initial values of design variables (increased thickness)
volfrac=0.3; % volumn fraction of solid material (upper bound)
E=2e11;nu=0.3; rho=7800; % Young's modulus, Poisson's ratio, density
mass=1e5; % concentrated mass
n_ord=1;n_ord_all=6; % n-th vibration frequency, all calculated frequencies
dgt0=5; % significant digit of sens.
scl_old=100; % initial set of scale factor for obj.
p=6;lmd=100; % power of super ellipsoid, power of KS aggregation
iter=1;maxiter=500; % initial and maximum number of iterations
objVr5=1.0; % initial relative variat. of obj. last 5 iterations
alpha=1e-6; mid_num=30; rank=0.5; ini=0.02; ep_old=0.3; % Parameter settings
% for the epsilon variation criterion
loop=[];var_obj=[];
% SEC 2): SETTING OF FE DISCRETIZATION
nEle=nelx*nely; % number of finite elements
nNod=(nelx+1)*(nely+1); % number of nodes
nDof=2*(nelx+1)*(nely+1); % number of degree of
% freedoms
EW=DW/nelx;EH=DH/nely; % length and width of
% finite elements
minSz=min([EW,EH]); % minimum size of finite
% elements
Ke=Ke_tril(E,nu,EW,EH,h); % upper triangular of ele.
% stiffness
KE(tril(ones(8))==1)=Ke';
KE=reshape(KE,8,8); KE=KE+KE'-diag(diag(KE)); % full elemental stiffness
% matrix
Me=Me_tril(rho,EW,EH,h); % upper triangular of ele.
% mass
ME(tril(ones(8))==1)=Me';
ME=reshape(ME,8,8);ME=ME+ME'-diag(diag(ME)); % full elemental mass
% matrix
nodMat=int32(reshape(1:nNod,1+nely,1+nelx)); % matrix of nodes numbers
% (int32)
edofVec=reshape(2*nodMat(1:end-1,1:end-1)-1,nEle,1);
edofMat=edofVec+int32([0 1 2*nely+[2 3 4 5] 2 3]); % connectivity matrix
eleNodesID=edofMat(:,2:2:8)./2;
[sI,sII]=deal([]);
for j=1:8
    sI=cat(2,sI,j:8);
    sII=cat(2,sII,repmat(j,1,8-j+1));
end
[iK,jK]=deal(edofMat(:,sI)',edofMat(:,sII)');
Iar0=sort([iK(:),jK(:)],2,'descend'); clear iK jK; % reduced assembly indexing
[x,y]=meshgrid(EW*(-nelx/2:nelx/2),EH*(-nely/2:nely/2)); % coordinates of
% nodal points
LSgrid.x=x(:);LSgrid.y=y(:);
volNod=sparse(double(eleNodesID(:)),1,1/4); % weight of each node in volume
% calculation
% SEC 3): LOADS, DISPLACEMENT BOUNDARY CONDITIONS
fixNd=union(1:nely+1,nelx*(nely+1)+1:(nelx+1)*(nely+1));
fixDof=[2*fixNd-1,2*fixNd];
fixEle=union(1:nely,(nelx-1)*nely+1:nely*nelx);
freeDof=setdiff(1:nDof,fixDof);

massNd=(nely+1)*nelx/2+nely/2+1;massDof=2*massNd-1:2*massNd;

massEle=[(nelx-1)*nely/2 (nelx-1)*nely/2+1 (nelx+1)*nely/2 (nelx+1)*nely/2+1];

% SEC 4): INITIAL SETTING OF COMPONENTS

xO=xInt:2*xInt:DW; yO=yInt:2*yInt:DH; % coordinates of initial components'
% center

xn=length(xO);yn=length(yO); % num. of component along x

xO=kron(xO,ones(1,2*yn));yO=repmat(kron(yO,ones(1,2)),1,xn); % full
% coordinates vector

N=length(xO); % total number of components

L=repmat(vInt(1),1,N); % vector of half-length

t1=repmat(vInt(2),1,N); % vector of 1st half-width

t2=repmat(vInt(3),1,N); % vector of 2nd half-width

theta=repmat([vInt(4) -vInt(4)],1,N/2); % vector of inclined angle

dd=[xO-EW*nelx/2; yO-EH*nely/2;L;t1;t2;theta]; % design variable vector

nDsvb=length(dd(:)); % number of all design variables

nEhcp=nDsvb/N; % number design variables each component

actComp=1:N; % initial set of active components

actDsvb=1:nDsvb; % initial set of active design variables

nNd=0;PhiNd=[]; % number of non-design patch and its TDF matrix

allPhi=[zeros(nNod,N) PhiNd]; % initialized TDF matrix

% SEC 5): SETTING OF MMA

m=1;c=1000*ones(m,1);d=zeros(m,1); % MMA parameters

a0=1;a=zeros(m,1);

xval=dd(:);xold1=xval;xold2=xval;

xmin=[-DW/2;-DH/2;minSz;minSz;minSz;-pi]; % lower bounds

xmax=[DW/2;DH/2;sqrt(DW^2+DH^2)/2;0.5*min(DW,DH)*[1;1];pi]; % upper bounds

xmin=repmat(xmin,N,1);xmax=repmat(xmax,N,1);

low=xmin; upp=xmax;

% Move limits for MMA
move = 0.1; % maximum move per iteration (relative to bounds)

% Store best valid design for recovery
xval_best = xval;
Lamb1_best = 0;

% SEC 6): OPTIMIZATION LOOP

OBJ=zeros(1,maxiter);
CONS=zeros(1,maxiter);

while objVr5>1e-4 && iter<=maxiter

    epsilon=ini+ep_old-(1/ep_old+exp(-rank*(iter-mid_num)))^-1;

    % LP 1): Generating TDFs and their derivatives

    allPhiDrv=sparse(nNod,nDsvb);

    for i=actComp % calculating TDF of the active MMCs

        [allPhi,allPhiDrv,xval,actComp,actDsvb]=...
            calc_Phi(allPhi,allPhiDrv,xval,i,LSgrid,p,nEhcp,epsilon,actComp,actDsvb,minSz);

    end

    allPhiAct=[allPhi(:,actComp) PhiNd]; % TDF matrix of active components

    temp=exp(lmd*allPhiAct);

    Phimax=max(-1e3,log(sum(temp,2))/lmd); % global TDF using K-S aggregation

    allPhiDrvAct=allPhiDrv(:,actDsvb);

    Phimaxdphi=kron(temp(:,1:length(actComp))./(sum(temp,2)+eps),ones(1,nEhcp));

    PhimaxDrvAct=Phimaxdphi.*allPhiDrvAct; % nodal sensitivity of global TDF

    % LP 2): Plotting current design
    figure(1);
    clf;
    contourf(reshape(x,[nely+1,nelx+1]),reshape(y,[nely+1,nelx+1]),...
        reshape(Phimax,[nely+1,nelx+1]),[0,0],'LineColor','none');
    colormap([1 1 1; 0 0 0]); % white background, black material
    axis equal; axis([-DW/2 DW/2 -DH/2 DH/2]);
    title(sprintf('Iteration %d',iter));
    drawnow;

    % LP 3): Finite element analysis

    H=Heaviside(Phimax,alpha,epsilon);

    den=sum(H(eleNodesID),2)/4;

    Psi1=zeros(nDof,1);

    bw_den=reshape(den,nely,nelx)>alpha;
    bw_den_sol=reshape(bwlabel(bw_den,4),nelx*nely,1);

    bw_fix=unique(bw_den_sol(fixEle));
    bw_mass=nonzeros(unique(bw_den_sol(massEle)));

    clear bw_den

    % Check if any mass region is connected to a fixed region
    hasLoadPath = false;
    connectedLabel = [];
    if ~isempty(bw_mass)
        for ilab = 1:length(bw_mass)
            if ismember(bw_mass(ilab), bw_fix)
                hasLoadPath = true;
                connectedLabel = bw_mass(ilab);
                break;
            end
        end
    end

    if hasLoadPath

        struct=1;

        eleLft=find(bw_den_sol==connectedLabel);

        denSld=zeros(size(den));denSld(eleLft)=den(eleLft);

        edofMatLft=edofMat(eleLft,:);

        freedofLft=setdiff(edofMatLft,fixDof); % retained DOFs for FEA

        [iK1,jK1]=deal(edofMatLft(:,sI)',edofMatLft(:,sII)');

        Iar=sort([iK1(:),jK1(:)],2,'descend'); clear iK1 jK1; % new reduced
        % assembly indexing

        sK=reshape(Ke(:)*denSld(eleLft)',length(Ke)*length(eleLft),1);
        K=fsparse(Iar(:,1),Iar(:,2),sK,[nDof,nDof]);K=K+K'-diag(diag(K));
        K=K+fsparse(1:nDof,1:nDof,eps*ones(1,nDof),[nDof,nDof]); %
        % regularization of disconnected component
        sM=reshape(Me(:)*denSld(eleLft)',length(Me)*length(eleLft),1);
        M=fsparse(Iar(:,1),Iar(:,2),sM,[nDof,nDof]);M=M+M'-diag(diag(M));
        M(massDof,massDof)=M(massDof,massDof)+mass*eye(length(massDof)); % add
        % concentrated mass
        M=M+fsparse(1:nDof,1:nDof,eps*ones(1,nDof),[nDof,nDof]); %
        % regularization of disconnected component
        [PsiAll,LambAll]=eigs(K(freedofLft,freedofLft),M(freedofLft,freedofLft),n_ord_all,'sm');
        Lamb1=LambAll(n_ord,n_ord);
        Psi1(freedofLft)=PsiAll(:,n_ord);Psi1=Psi1/sqrt(Psi1.'*M*Psi1);

        % Update best design if this is better
        if Lamb1 > Lamb1_best
            Lamb1_best = Lamb1;
            xval_best = xval;
        end
    else % no load path - use full domain with penalty
        disp('WARNING!!! NO loading path is founded - applying penalty!!!');
        struct=0;

        % Use full domain FEA
        sK=reshape(Ke(:)*den(:)',length(Ke)*nEle,1);
        K=fsparse(Iar0(:,1),Iar0(:,2),sK,[nDof,nDof]);K=K+K'-diag(diag(K));
        K=K+fsparse(1:nDof,1:nDof,eps*ones(1,nDof),[nDof,nDof]); % regularization
        sM=reshape(Me(:)*den(:)',length(Me)*nEle,1);
        M=fsparse(Iar0(:,1),Iar0(:,2),sM,[nDof,nDof]);M=M+M'-diag(diag(M));
        M(massDof,massDof)=M(massDof,massDof)+mass*eye(length(massDof)); % add
        % concentrated mass
        M=M+fsparse(1:nDof,1:nDof,eps*ones(1,nDof),[nDof,nDof]); % regularization
        [PsiAll,LambAll]=eigs(K(freeDof,freeDof),M(freeDof,freeDof),n_ord_all,'sm');
        Lamb1=LambAll(n_ord,n_ord);
        Psi1(freeDof)=PsiAll(:,n_ord);Psi1=Psi1/sqrt(Psi1.'*M*Psi1);

        % Apply heavy penalty for disconnection (use fraction of best value)
        if Lamb1_best > 0
            Lamb1 = 0.1 * Lamb1_best; % Penalize disconnected state
        end
    end
    OBJ(iter)=Lamb1; % scaled objective function
    fval=sum(den)*EW*EH/(DW*DH)-volfrac; CONS(iter) = fval + volfrac; % volume
    % constraint
    % LP 4): Sensitivity analysis
    df0dx=zeros(1,nDsvb);dfdx=zeros(1,nDsvb);
    delta_H=3*(1-alpha)/(4*epsilon)*(1-Phimax.^2/(epsilon^2));
    delta_H(abs(Phimax)>epsilon)=0; % derivative of nodal density to nodal TDF
    if struct == 1
        [a1_w,a2_w]=calc_f0val(K,M,freedofLft,Psi1,nDof,Lamb1);
    else
        [a1_w,a2_w]=calc_f0val(K,M,freeDof,Psi1,nDof,Lamb1);
    end
    sens=sum((a1_w(edofMat)*(KE-Lamb1.*ME)).*Psi1(edofMat)-0.5*(a2_w*Psi1(edofMat)*ME).*Psi1(edofMat),2);
    sEner=sens*ones(1,4)/4; sensNod=sparse(double(eleNodesID(:)),1,sEner(:));
    % nodal form
    df0dx(actDsvb)=(sensNod.*delta_H)'*PhimaxDrvAct; % sensitivity of
    % objective function
    dfdx(actDsvb)=(volNod.*delta_H)'*PhimaxDrvAct*EW*EH/(DW*DH); % sensitivity
    % of volume constraint
    dgt=dgt0-floor(log10([max(abs(df0dx(:)))+eps max(abs(dfdx(:)))+eps])); %
    % significant digits for sens. truncation
    scl=max(abs(df0dx));
    if scl < eps || ~isfinite(scl)
        scl = scl_old; % Use previous scale if current is invalid
    else
        scl_old = scl;
    end
    f0val=-Lamb1/scl;
    df0dx=round(df0dx*10^dgt(1))/10^dgt(1)/scl; % truncated scaled objective
    % sensitivity
    dfdx=round(dfdx*10^dgt(2))/10^dgt(2); % truncated constraint sensitivity
    % LP 5): Updating design variables
    % Apply move limits
    xmin_move = max(xmin, xval - move*(xmax-xmin));
    xmax_move = min(xmax, xval + move*(xmax-xmin));

    [xmma,~,~,~,~,~,~,~,~,low,upp]=mmasub(m,nDsvb,iter,xval,...
        xmin_move,xmax_move,xold1,xold2,f0val,df0dx(:),fval,dfdx,low,upp,a0,a,c,d);
    xold2=xold1;xold1=xval;xval=xmma; % design variable's update
    if iter>=5 && fval/volfrac<1e-4
        objVr5=abs(max(abs(OBJ(iter-4:iter)-mean(OBJ(iter-4:iter))))/mean(OBJ(iter-4:iter)));
    end

    disp(['It.: ' sprintf('%4i\t',iter) 'Obj.: ' sprintf('%6.3f\t',Lamb1) 'Vol.: ' ...
        sprintf('%6.4f\t',fval) 'Ch.: ' sprintf('%6.4f\t',objVr5)]);
    % figure(2);
    var_obj=[var_obj Lamb1];loop=[loop iter];
    % yyaxis left
    % plot(loop,var_obj,'DisplayName','w1');
    % xlabel('Step');ylabel('Frequency');
    % yyaxis right
    % plot(loop,CONS,'DisplayName','vol')
    % ylabel('Vol');ylim([0 1]);
    % legend('show','location','southeast');
    % drawnow
    iter=iter+1;
end

function [allPhi,allPhidrv,xval,actComp,actDsvb] = ...
    calc_Phi(allPhi,allPhidrv,xval,i,LSgrid,p,nEhcp,epsilon,actComp,actDsvb,minSz)
dd = xval((i-1)*nEhcp+1:i*nEhcp);
x0 = dd(1); y0 = dd(2); L = dd(3) + eps; t1 = dd(4); t2 = dd(5);
st = sin(dd(6)); ct = cos(dd(6));
x1 = ct*(LSgrid.x-x0) + st*(LSgrid.y-y0) + eps; % local x of
% each grid
y1 = -st*(LSgrid.x-x0) + ct*(LSgrid.y-y0) + eps; % local y of
% each grid
l = (t1+t2)/2 + (t2-t1)/2/L*x1 + eps;
temp = abs(x1).^p/L^p + abs(y1).^p./l.^p;
allPhi(:,i) = 1 - temp.^(1/p); % TDF of i-th
% component
if dd(4)/minSz < 1.01 && dd(5)/minSz < 1.01
    % disp(['The ' sprintf('%i',i) '-th component is too small! DELETE it!!!']);
    allPhi(:,i) = -1e3;
    xval((i-1)*nEhcp+(4:5)) = 0;
    actComp = setdiff(actComp,i);
    actDsvb = setdiff(actDsvb,nEhcp*i-nEhcp+1:nEhcp*i);
elseif min(abs(allPhi(:,i))) >= epsilon
    % disp(['The ' sprintf('%i',i) '-th component is too small! DELETE it!!!']);
    allPhi(:,i) = -1e3;
    xval((i-1)*nEhcp+(4:5)) = 0;
    actComp = setdiff(actComp,i);
    actDsvb = setdiff(actDsvb,nEhcp*i-nEhcp+1:nEhcp*i);
else % calculate derivatives
    dx1 = [-ct+0.0*x1 -st+0.0*x1 0.0*x1 0.0*x1 0.0*x1 y1]; % variation
    % of x'
    dy1 = [st+0.0*y1 -ct+0.0*y1 0.0*y1 0.0*y1 0.0*y1 -x1]; % variation
    % of y'
    dldx1 = (t2-t1)/2/L;
    dldv1 = 0.0*l; dldv2 = 0.0*l; dldv6 = 0.0*l; % dldx0,
    % dldy0, dldtheta
    dldv3 = -(t2-t1)/2*x1/L^2; dldv4 = 1/2 - x1/2/L; dldv5 = 1/2 + x1/2/L;
    dl = [dldv1 dldv2 dldv3 dldv4 dldv5 dldv6] + repmat(dldx1,1,nEhcp).*dx1; %
    % variation of width
    dpdx1 = -(temp).^(1/p-1).*((x1/L).^(p-1)/L); % dphi/dx

    dpdy1 = -(temp).^(1/p-1).*((y1./l).^(p-1)./l); % dphi/dy

    dpdL = zeros(size(dx1)); dpdL(:,3) = (temp).^(1/p-1).*(x1/L).^p/L; % [0 0
    % dphi/dL 0 0 0]
    dpdl = (temp).^(1/p-1).*((y1./l).^p./l);
    allPhidrv(:,nEhcp*(i-1)+1:nEhcp*i) = repmat(dpdx1,1,nEhcp).*dx1 +...
        repmat(dpdy1,1,nEhcp).*dy1 + dpdL + repmat(dpdl,1,nEhcp).*dl;
end
end

function [a1,a2]=calc_f0val(K,M,freeDof,Psi,nDof,Lamb)

L11=K(freeDof,freeDof)-Lamb*M(freeDof,freeDof);

L12=-M(freeDof,freeDof)*Psi(freeDof);L21=L12';

LL = [L11,L12;
      L21,0];

vec_a=LL\[zeros(length(freeDof),1);1];

a1=zeros(nDof,1);a1(freeDof)=vec_a(1:end-1);a2=vec_a(end);

end

function [H] = Heaviside(phi,alpha,epsilon)

H = 3*(1-alpha)/4*(phi/epsilon-phi.^3/(3*(epsilon)^3)) + (1+alpha)/2;

H(phi>epsilon) = 1;

H(phi<-epsilon) = alpha;

end


function Ke=Ke_tril(E,nu,a,b,h)

k1=[-1/6/a/b*(nu*a^2-2*b^2-a^2), 1/8*nu+1/8,-1/12/a/b*(nu*a^2+4*b^2-a^2), 3/8*nu-1/8,...
    1/12/a/b*(nu*a^2-2*b^2-a^2),-1/8*nu-1/8, 1/6/a/b*(nu*a^2+b^2-a^2), -3/8*nu+1/8];

k2=[-1/6/a/b*(nu*b^2-2*a^2-b^2), 1/8*nu+1/8,-1/12/a/b*(nu*b^2+4*a^2-b^2), 3/8*nu-1/8,...
    1/12/a/b*(nu*b^2-2*a^2-b^2),-1/8*nu-1/8, 1/6/a/b*(nu*b^2+a^2-b^2), -3/8*nu+1/8];

Ke=E*h/(1-nu^2)*...
    [k1(1);k1(2);k1(3);k1(4);k1(5);k1(6);...
     k1(7);k1(8);k2(1);k2(8);k2(7);k2(6);...
     k2(5);k2(4);k2(3);k1(1);k1(6);k1(7);...
     k1(4);k1(5);k1(2);k2(1);k2(8);k2(3);...
     k2(2);k2(5);k1(1);k1(2);k1(3);k1(4);...
     k2(1);k2(8);k2(7);k1(1);k1(6);k2(1)];

end


function Me=Me_tril(rho,EW,EH,h)

A=EW*EH;

Me=(rho*A*h/36)*...
    [4;0;2;0;1;0;...
     2;0;4;0;2;0;...
     1;0;2;4;0;2;...
     0;1;0;4;0;2;...
     0;1;4;0;2;0;...
     4;0;2;4;0;4];

end
