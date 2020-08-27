%SLOW pruning algorithm
function spec_mat=expmeldor_exchange_prune(cA,ksym,Rpll,B0,ntmix,angle1,angle2)
PRUNE_TOL=1e-2;

rngeMHz=333;%500;%omarr goes from -rngeMHz to +rngeMHz
NGRD=128;%number of grid points in f1,f2

%fixed for the paper
gxx=2.0087;
gyy=2.0057;
gzz=2.0021;
shiftroffdiag=2.0;%G %1.0 %shiftroffdiag added to off-diag space only

tmix=ntmix*1e-3;%take to microsecs
cfact=1e-6*mean([gxx,gyy,gzz])*9.2731e-21/1.05443e-27; %changes with frequency(?)
strng=strcat(num2str(B0),'_',num2str(Rpll));    

%define kf and kr
ksymG=ksym/(cfact/(2*pi));%per microsecond to Gauss
k_tot=ksymG/sqrt(cA*(1-cA));
kf=(1-cA)*k_tot;kr=cA*k_tot;

%define stvx
stvxA=load(strcat('stvech_',strng,'_',angle1,'.stvx'));
stvxB=load(strcat('stvech_',strng,'_',angle2,'.stvx'));
stvx=[sqrt(cA)*(stvxA(:,1)+1i*stvxA(:,2));sqrt(1-cA)*(stvxB(:,1)+1i*stvxB(:,2))];

omarr=2*pi*linspace(-rngeMHz,rngeMHz,NGRD);
omarrG=omarr/cfact; %calculate in Gauss units

%define off-diag space matrix
a=load(strcat('matlab_realh_',strng,'_',angle1,'.mtxx'));
matxrA=spconvert(a);
a=load(strcat('matlab_realh_',strng,'_',angle2,'.mtxx'));
matxrB=spconvert(a);
a=load(strcat('matlab_imagh_',strng,'_',angle1,'.mtxx'));
matxiA=spconvert(a);
a=load(strcat('matlab_imagh_',strng,'_',angle2,'.mtxx'));
matxiB=spconvert(a);

%calculate shape of matrix, define the matrix (should have kf,kr in it)
ndimo=size(matxrA,1);
ndimo_orig=ndimo;
size(matxrA)
size(matxrB)

matx=[matxrA+1i*matxiA-kf*speye(ndimo_orig),ksym*speye((ndimo_orig));ksym*speye(ndimo_orig),matxrB+1i*matxiB-kr*speye(ndimo_orig)]+shiftroffdiag*speye(2*ndimo_orig);
%issparse(matx)
ndimo=2*ndimo;
nnz(matx)

disp('Issparse');
disp(issparse(matx));

%define diag space matrix
a=load(strcat('matlab_realh_',strng,'_',angle1,'.mtxz'));
matzrA=spconvert(a);
a=load(strcat('matlab_realh_',strng,'_',angle2,'.mtxz'));
matzrB=spconvert(a);
a=load(strcat('matlab_imagh_',strng,'_',angle1,'.mtxz'));
matziA=spconvert(a); 
a=load(strcat('matlab_imagh_',strng,'_',angle2,'.mtxz'));
matziB=spconvert(a); 

%calculate shape of matrix, define the matrix (should have kf,kr in it)
ndimd_orig=size(matzrA,1);
matz=[matzrA+1i*matziA-kf*speye((ndimd_orig)),ksym*speye((ndimd_orig));ksym*speye((ndimd_orig)),matzrB+1i*matziB-kr*speye((ndimd_orig))];


r=0*omarrG;%array of rms errors in resolvents

tic
stvzchk=zeros(ndimo,1);
%should I make this loop parfor?
preprunecwspec = [];
for ii=0:int32(numel(omarrG)/5-1)
    j=5*ii+1;%pick ~nged/5 points
    invec = (matx+1i*speye(ndimo)*omarrG(j))\stvx;
    gmres(matx+1i*speye(ndimo)*omarrG(j))\stvx;
    res = invec/(stvx'*invec);
    stvzchk = stvzchk + (abs(res) > stvzchk).*(abs(res)-stvzchk); %note that we are summing over all the abs values, no mean/rms here!
    preprunecwspec = [preprunecwspec;stvx'*invec];
end
toc
PINDX=find(stvzchk>=PRUNE_TOL); %PINDX has all the pruned basis vectors in the off-diag space
disp('reduction');
disp(length(PINDX)/ndimo)
disp('max(PINDX)/ndimo');
disp(max(PINDX)/ndimo);
disp([cA,ksym,Rpll,B0,ntmix,angle1,angle2, max(PINDX)/ndimo])
%prune matx
matx=matx(PINDX,PINDX);%PRUNE THE OFF_DIAG MATRIX
ndimo=size(matx,1);%redefine ndimo, orig ndimo in ndimo_orig
stvx=stvx(PINDX);

stvznop=zeros(ndimo,numel(omarrG));

tic

parfor i=1:numel(omarrG)
    stvznop(:,i)=(matx+1i*speye(ndimo)*omarrG(i))\stvx;
    r(i)=rms((matx+1i*speye(ndimo)*omarrG(i))*stvznop(:,i)-stvx);
end
rms(r)
toc

postprunecwspec=[];
for ii=0:int32(numel(omarrG)/5-1)
    j=5*ii+1;%pick the same ~ngrd/5 points
    postprunecwspec=[postprunecwspec;stvx'*stvznop(:,j)];
end

disp('PRUNE SUCCESS');
disp(rms((postprunecwspec-preprunecwspec)./preprunecwspec))
if(rms((postprunecwspec-preprunecwspec)./preprunecwspec)>0.01)
    disp('PRUNE BAD!!!!!!!');
end
tic
pprop=load(strcat('pproph_',strng,'_Exch.txt'));
pid=load(strcat('pidh_',strng,'_Exch.txt'));
indx=[];
for k=1:ndimo_orig
    if(pid(k)==1)
        indx=[indx;k];
    else
        if(pid(k)==2)
            indx=[indx;k;k];
        end
    end
end

pp=spconvert([(1:ndimd_orig)',indx,pprop]);
pp=blkdiag(pp,pp);%pulse propagators, same for both components A and B

%simple way to find out which of the diag space basis vectors matter
test=zeros(2*ndimo_orig,1);
test(PINDX)=1;
testz=pp*test;
PINDXZ=find(testz~=0);

%prune matxz
matz=matz(PINDXZ,PINDXZ);

stvz=pp(PINDXZ,PINDX)*stvznop;%tranform off-diag to diag
%spec_mat=stvz'*expm(-cfact*tmix*matz)*stvz;%calculate the 2d ELDOR spectrum
spec_mat=zeros(NGRD,NGRD);
parfor i=1:NGRD
    spec_mat(:,i)=stvz'*expv(-cfact,tmix*matz,stvz(:,i));%calculate the 2d ELDOR spectrum
end

dlmwrite(strcat('Paper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',num2str(angle2),'deg_ksym',num2str(ksym*1e6),'_r.txt'),real(spec_mat));
dlmwrite(strcat('Paper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',num2str(angle2),'deg_ksym',num2str(ksym*1e6),'_i.txt'),imag(spec_mat));

%plot!!
% figure;
% contour(-omarr/(2*pi),-omarr/(2*pi),abs(spec_mat)/max(max(abs(spec_mat))),30);%have to invert omarr
% xlabel('f_1 (MHz)');
% ylabel('f_2 (MHz)');
% set(gca,'FontSize',16);
% title(strcat('k_{sym}=',num2str(ksym*1e6),' s^{-1}'),'Interpreter','latex','FontSize',20);
figure;
%mesh(-omarr,-omarr,abs(spec_mat)/max(max(abs(spec_mat))));%have to invert omarr
mesh(-omarr/(2*pi),-omarr/(2*pi),abs(spec_mat));%have to invert omarr
view([1,-1,0.5])
xlabel('f_1 (MHz)');
ylabel('f_2 (MHz)');
set(gca,'FontSize',16);
title(strcat('$k_{sym}=',num2str(ksym*1e6),'\ s^{-1}, T_{mix}=',num2str(ntmix),' ns, \beta_d=',angle2,'$'),'Interpreter','latex','FontSize',20);
saveas(gcf,strcat('Paper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',num2str(angle2),'deg_ksym',num2str(ksym*1e6),'.fig'));
saveas(gcf,strcat('Paper1exch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_betadjump',num2str(angle2),'deg_ksym',num2str(ksym*1e6),'.pdf'));
close(gcf);
%saveas(gcf,strcat('ksym',num2str(ksym*1e6),'_Tmix',num2str(ntmix),'_Betad',angle2,'.fig'));
%title({'new |S_{c-}|', strcat('T_{mix}=0 ns, R_{pll}=10^, inter-conformational angle',angle1,num2str(Rpll),' s^{-1}')},'Fontsize',16);
%dlmwrite(strcat('resolvents_Tmix',num2str(ntmix),'ns_',strng,'.out'),abs(spec_mat));
toc

end

