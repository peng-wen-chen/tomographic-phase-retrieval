function RHS=CG_RHS(Fh1,Fh2,Fh3,beta,input)
% output is the right-hand side of the CG normal equation
% P_{Z_n} FFT(K^* Fh)
rho=input.rho;
n=input.n;
s=input.s;
p=input.p;
p1=input.p1;
i_index=[-p1/2:p1/2]'+s;


beta1=beta.beta1;
beta2=beta.beta2;
beta3=beta.beta3;

for kk=1:p %  construct K^* Fh1, K^* Fh2, K^* Fh3, 
    for jj=1:p
    Kh1(1:p,jj,kk)=(exp(-i*2*pi/p*i_index*(beta1(:,1)*rho(jj)+beta1(:,2)*rho(kk))')*squeeze(Fh1(:,jj,kk)));% p*p
    Kh2(jj,1:p,kk)=(exp(-i*2*pi/p*i_index*(beta2(:,1)*rho(jj)+beta2(:,2)*rho(kk))')*transpose(Fh2(jj,:,kk)));% p*p
    Kh3(jj,kk,1:p)=(exp(-i*2*pi/p*i_index*(beta3(:,1)*rho(jj)+beta3(:,2)*rho(kk))')*squeeze(Fh3(jj,kk,:)));% p*p
    end
end  

for ii=1:p
FKh1(ii,:,:)=fftshift(ifft2(ifftshift(squeeze(Kh1(ii,:,:)))))*p;
FKh2(:,ii,:)=fftshift(ifft2(ifftshift(squeeze(Kh2(:,ii,:)))))*p;
FKh3(:,:,ii)=fftshift(ifft2(ifftshift(squeeze(Kh3(:,:,ii)))))*p;
end

FKh_sum=FKh1+FKh2+FKh3;
% from here to compute the gradient projection
% f part
 clear Fg
 % Apply P_{Z_n}, the support constraint projection
 RHS=zeros(p,p,p);
 ind1=input.ind1;
ind2=input.ind2;
 RHS(ind1:ind2,ind1:ind2,ind1:ind2)=FKh_sum(ind1:ind2,ind1:ind2,ind1:ind2);

end
