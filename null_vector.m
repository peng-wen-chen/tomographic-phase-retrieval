function [fs, cor] = null_vector(weight,fs_ini, input)
mu=input.mu;
n=input.n;p0=input.p0;
f0=input.fs;
m=input.m;

ind1=weight.data1;
ind2=weight.data2;
ind3=weight.data3;

norm_f0=norm(f0(:));
IterMax=input.Nuliter;
 fs=fs_ini;
%%
for iter=1:IterMax % power iteration
[h1,h2,h3]=Rmap(fs, input);
gF1=zeros(m,p0,p0);
gF2=zeros(p0,m,p0);
gF3=zeros(p0,p0,m);

for ii=1:m
    gF1(ii,:,:)=Qmap(h1(ii,:,:), squeeze(mu(ii,:,:)),p0);
    gF2(:,ii,:)=Qmap(h2(:,ii,:), squeeze(mu(ii,:,:)),p0);
    gF3(:,:,ii)=Qmap(h3(:,:,ii), squeeze(mu(ii,:,:)),p0);
end
  z1=gF1;   z2=gF2;    z3=gF3;

  
    index1=input.ind1;
    index2=input.ind2;
f0out=threeD_to_twoD(fs(index1:index2,index1:index2,index1:index2),input);
figure(1); clf; imshow(abs(f0out)); % show the moduli of null vector
cor(iter)=abs(fs(:)'*f0(:))/norm(fs(:))/norm(f0(:));
sprintf('At %d th iteration, correlation is %15.4f',iter,cor(iter))

    bu=ind1.*z1; 
for ii=1:m
        h1(ii,:,:)=Qmap_inv(bu(ii,:,:), squeeze(mu(ii,:,:)),input );
end

bu=ind2.*z2;
for ii=1:m
   h2(:,ii,:)=Qmap_inv(bu(:,ii,:), squeeze(mu(ii,:,:)),input);
end
bu=ind3.*z3;
for ii=1:m
     h3(:,:,ii)=Qmap_inv(bu(:,:,ii), squeeze(mu(ii,:,:)),input);
end
[fs]=Rmap_inv(h1, h2,h3,fs, input);


 
 fs = P_omega(fs, input); % if phantom, projection on Q1 is used

fs=fs/norm(fs(:))*norm_f0; %, norm(fs(:)));
end
end

