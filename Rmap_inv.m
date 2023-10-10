function [fs]=Rmap_inv(Fh1,Fh2,Fh3,f_ini, input)
p=input.p;
n=input.n;
s=input.s;
beta.beta1=input.beta1;
beta.beta2=input.beta2;
beta.beta3=input.beta3;


m=input.m;

for ii=1:m
h1(ii,:,:)=fftshift(fft2(ifftshift(squeeze(Fh1(ii,:,:)))))/p;
h2(:,ii,:)=fftshift(fft2(ifftshift(squeeze(Fh2(:,ii,:)))))/p;
h3(:,:,ii)=fftshift(fft2(ifftshift(squeeze(Fh3(:,:,ii)))))/p;


end


A_h=CG_RHS(h1,h2,h3,beta,input); % A_h is the RHS of CG normal equation


 fs = SolveNormalCG(A_h , f_ini,input);%p,T_hat,n);


end
