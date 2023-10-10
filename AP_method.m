function [fs, err]=AP_method(b,fs_ini,input)
% The code is AP method, which reconstructs fs from b, starting with
% initialization fs_ini

beta.beta1=input.beta1;
beta.beta2=input.beta2;
beta.beta3=input.beta3;
mu=input.mu;
n=input.n;
p0=input.p0;
p=input.p;
m=input.m;
p1=input.p1;

b1=b.data1;
 b2=b.data2;
 b3=b.data3;

IterMax=input.APiter;
fs = P_omega(fs_ini, input);

%% construct u from fs

[h1,h2,h3]=Rmap(fs, input);

for ii=1:m
gF1(ii,:,:)=Qmap(h1(ii,:,:), squeeze(mu(ii,:,:)),p0);
    gF2(:,ii,:)=Qmap(h2(:,ii,:), squeeze(mu(ii,:,:)),p0);
    gF3(:,:,ii)=Qmap(h3(:,:,ii), squeeze(mu(ii,:,:)),p0);
end

 u1=gF1./(abs(gF1)+1e-6);
 u2=gF2./(abs(gF2)+1e-6);
 u3=gF3./(abs(gF3)+1e-6);



%% AP algorithm
for iter=1:IterMax
    % Qmap_inv
bu=b1.*u1;
for ii=1:m
    h1(ii,:,:)=Qmap_inv(bu(ii,:,:), squeeze(mu(ii,:,:)),input);
end
 bu=b2.*u2;
for ii=1:m
        h2(:,ii,:)=Qmap_inv(bu(:,ii,:), squeeze(mu(ii,:,:)),input);
end
 bu=b3.*u3;
for ii=1:m
        h3(:,:,ii)=Qmap_inv(bu(:,:,ii), squeeze(mu(ii,:,:)),input);
end
% Rmap_inv
[fs]=Rmap_inv(h1,h2,h3,fs, input);
fs = P_omega(fs, input);
[h1,h2,h3]=Rmap(fs, input);
for ii=1:m
gF1(ii,:,:)=Qmap(h1(ii,:,:), squeeze(mu(ii,:,:)),p0);
     gF2(:,ii,:)=Qmap(h2(:,ii,:), squeeze(mu(ii,:,:)),p0);
    gF3(:,:,ii)=Qmap(h3(:,:,ii), squeeze(mu(ii,:,:)),p0);
end

 u1=gF1./(abs(gF1)+1e-6);
 u2=gF2./(abs(gF2)+1e-6);
 u3=gF3./(abs(gF3)+1e-6);

 ind1=input.ind1;
 ind2=input.ind2; 
f0out=threeD_to_twoD(fs(ind1:ind2,ind1:ind2,ind1:ind2),input);
figure(1)
imshow(abs(f0out))
err(iter)=norm([abs(gF1(:))-b1(:);...
    abs(gF2(:))-b2(:); ...
    abs(gF3(:))-b3(:)])/norm([b1(:);...
    b2(:);b3(:)]);
disp('[iter, err]=')
[iter err(iter)]

end
end