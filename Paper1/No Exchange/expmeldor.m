%Paper check:
%check the Shell script file 
%check the rr.sh driver file for the Shell script
%What am I doing new here? Matrices OK by default, pp I made a correction,
%time-freq-mag field conversion: need to confirm
%Check the flow, whether MATLAB equations are correct (the flowchart of the "new algo" in the paper)
%In the pruning program check pruning

%fixed for the paper
gxx=2.0087;
gyy=2.0057;
gzz=2.0021;
shiftroffdiag=2.0;%G %1.0


%loop over mix time
for ntmix=[0]%[0,50,100] %[50,100,150,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000]%0.000001
%loop over B0 fields Gauss
for B0=34050 %3300%[3300,34050]
%loop over Rpll; assuming that Rprp is 1/2 Rpll
for Rpll=[5]%[5,6]%9%5 %[5.5 5] %[9 8.5 8 7.5 7 6.5 6]
%count time
tic
%128 grid points in f1, f2 both 128x128 2d grid
NGRD=128;

%convert nanosecs to microsecs
tmix=ntmix*1e-3;%take to microsecs

%take microsecs to 1/Gauss
cfact=1e-6*mean([gxx,gyy,gzz])*9.2731e-21/1.05443e-27;%changes with frequency(?), ans should be 17.6372 for the default gii vals

%file ID string, expect the same in the nlspmc files that we generate
strng=strcat(num2str(B0),'_',num2str(Rpll));

%starting vector file
stvx=load(strcat('stvech_',strng,'.stvx'));
stvx=stvx(:,1)+1i*stvx(:,2); %real and im parts combine

a=load(strcat('matlab_realh_',strng,'.mtxx'));
matxr=spconvert(a);

a=load(strcat('matlab_imagh_',strng,'.mtxx'));
matxi=spconvert(a);

a=load(strcat('matlab_realh_',strng,'.mtxz'));
matzr=spconvert(a);

a=load(strcat('matlab_imagh_',strng,'.mtxz'));
matzi=spconvert(a); 

ndimo=size(matxr,1);%off-diag space size
ndimd=size(matzr,1);%diag space size

%matrices defined!
matx=matxr+shiftroffdiag*sparse(eye(size(matxr,1)))+1i*matxi; %shiftroffdiag added in diag space
matz=matzr+1i*matzi; %no shift added in diag space

%-333 to 333 MHz plotting
rngeMHz=333;%500;
omarr=2*pi*linspace(-rngeMHz,rngeMHz,NGRD); %rad/s, it is going to be converted into Gauss units
omarrG=omarr/cfact;

%residual: just for diagnostics
r=0*omarrG;
%we multiply by the pp here: goal is to get identity matrix for pp'*pp

stvznop=zeros(ndimo,numel(omarrG));

%disp('Non-normality condition number:');
%disp(norm(full(matx*matx'-matx'*matx))/norm(full(matx))^2)
imagesc(log10(abs(full(matx))));colorbar;return;
parfor i=1:numel(omarrG)
    stvznop(:,i)=(matx+1i*speye(ndimo)*omarrG(i))\stvx; %resolvents
    r(i)=rms((matx+1i*speye(ndimo)*omarrG(i))*stvznop(:,i)-stvx);
end
rms(r)
toc

tic
%pulse propagator
pprop=load(strcat('pproph_',strng,'.txt'));
pid=load(strcat('pidh_',strng,'.txt'));
indx=[];
for k=1:ndimo
    if(pid(k)==1)
        indx=[indx;k];
    else
        if(pid(k)==2)
            indx=[indx;k;k];
        end
    end
end

pp=spconvert([(1:ndimd)',indx,pprop]); %creates an ndimdxndimo matrix
if(mean(abs(pp'*pp-2*eye(ndimo)))>1e-15)
    disp('PP ERROR!!!!!!!');
end
%propagate to diagonal space
stvz=pp*stvznop;

%calculate the spectrum, done!!!
parfor i=1:NGRD
spec_mat(:,i)=stvz'*expv(-cfact,tmix*matz,stvz(:,i));
end

%PLOT!!
%spec_mat=stvz'*expm(-cfact*tmix*matz)*stvz;
dlmwrite(strcat('Paper1noexch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_r.txt'),real(spec_mat));
dlmwrite(strcat('Paper1noexch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns_i.txt'),imag(spec_mat));
figure;
%contour(omarr,omarr,abs(spec_mat)/max(max(abs(spec_mat))),30);
mesh(-omarr/(2*pi),-omarr/(2*pi),abs(spec_mat)/sum(sum(abs(spec_mat))));
view([1,-1,0.5]);
set(gca,'Fontsize',16);
xlabel('f_1 (MHz)', 'Fontsize',16);
ylabel('f_2 (MHz)', 'Fontsize',16);
title({'new |S_{c-}|', strcat('T_{mix}=',num2str(ntmix),' ns, R_{pll}=10^',num2str(Rpll),' s^{-1}')},'Fontsize',16);
saveas(gcf,strcat('Paper1noexch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns.fig'));
saveas(gcf,strcat('Paper1noexch_spech_mat_',strng,'_Tmix',num2str(ntmix),'ns.pdf'));
close(gcf);
toc
end
end
end