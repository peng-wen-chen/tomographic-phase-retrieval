 function [h1,h2, h3]=Rmap(f, input)
% we update s from 0 to n/2
p=input.p;
n=input.n;
m=input.m;
s=input.s;
beta1=input.beta1;
beta2=input.beta2;
beta3=input.beta3;
% this is R from Omega to Z_p^2
% T1, T2, T3 is size of pb*p*p
% compute A from f(i,j) to h(beta, rho),h consist of T1 and T2
F1=zeros(p,p,p);

for ii=1:p
F1(ii,:,:)=fftshift(fft2(ifftshift(squeeze(f(ii,:,:)))))/p;
end
p1=input.p1;
% if rem(p,2)==1
%     p1=p-1;i_index=[-p1/2:p1/2]; % same as rho
% 
% else
%     p1=p;
%   
% end
  i_index=[-p1/2:p1/2]; % same as rho
for kk=1:p % rho line
    for jj=1:p
    Fh1(:,jj,kk)=(exp(i*2*pi/p*(beta1(:,1)*i_index(jj)+beta1(:,2)*i_index(kk))*(i_index+s))*squeeze(F1(:,jj,kk)));% p*p
    end
end  


for ii=1:m
h1(ii,:,:)=fftshift(ifft2(ifftshift(squeeze(Fh1(ii,:,:)))))*p;
end


F2=zeros(p,p,p);
for ii=1:p
F2(:,ii,:)=fftshift(fft2(ifftshift(squeeze(f(:,ii,:)))))/p;
end
for kk=1:p % rho line
    for jj=1:p
    Fh2(jj,:,kk)=(exp(i*2*pi/p*(beta2(:,1)*i_index(jj)+beta2(:,2)*i_index(kk))*(i_index+s))*transpose(squeeze(F2(jj,:,kk))));% p*p
    end
end 


for ii=1:m
h2(:,ii,:)=fftshift(ifft2(ifftshift(squeeze(Fh2(:,ii,:)))))*p;
end

F3=zeros(p,p,p);

for ii=1:p
F3(:,:,ii)=fftshift(fft2(ifftshift(squeeze(f(:,:,ii)))))/p;
end
for kk=1:p % rho line
    for jj=1:p
    Fh3(jj,kk,:)=(exp(i*2*pi/p*(beta3(:,1)*i_index(jj)+beta3(:,2)*i_index(kk))*(i_index+s))*(squeeze(F3(jj,kk,:))));
    end
end

for ii=1:m
h3(:,:,ii)=fftshift(ifft2(ifftshift(squeeze(Fh3(:,:,ii)))))*p;
end


end
