%check quality of eigenvectors
%zcondarr=zeros(6,1);
xcondarr=zeros(6,1);
cnt =0;
for Rpll=[5,6,7,8,9,10]
    matxr = spconvert(load(strcat('egalmatlab_realh_34050_',num2str(Rpll),'.mtxx')));
    matxi = spconvert(load(strcat('egalmatlab_imagh_34050_',num2str(Rpll),'.mtxx')));
    %matzr = spconvert(load(strcat('matlab_realh_34050_',num2str(Rpll),'.mtxz')));
    %matzi = spconvert(load(strcat('matlab_imagh_34050_',num2str(Rpll),'.mtxz')));
    matx = matxr + 1i* matxi;
    %matz = matzr + 1i* matzi;
    [evecx,~] = eig(full(matx));
    %[evecz,evalz] = eig(full(matz));
    cnt = cnt+1;
    %zcondarr(cnt)=norm(evecz)*norm(inv(evecz));
    xcondarr(cnt)=norm(evecx)*norm(inv(evecx));
    clear evecx;
    clear matx;
end
%figure;
%plot([5,6,7,8,9,10],zcondarr);
figure;
plot([5,6,7,8,9,10],xcondarr);
%figure;
%aa= abs(evecx.'*evecx-diag(diag(evecx.'*evecx)));
%inds=1:int32(size(evecx,1)/500):size(evecx,1);
%contour(aa(inds,inds))
%%
figure;
mesh(abs(evecz.'*evecz-diag(diag(evecz.'*evecz))));
figure;
mesh(abs(matx - evecx * evalx * evecx.'));


%%
stvx=[1;1;1;zeros(size(matx,1)-3,1)]/sqrt(3);
tic
x = gmres(matx+2*speye(size(matx,1)),stvx,5,1e-10,100);
%x = (matx+2*speye(size(matx,1)))\stvx;
%rms((matx+2*speye(size(matx,1)))*x-stvx)
toc