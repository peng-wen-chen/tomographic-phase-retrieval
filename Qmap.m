function gF=Qmap(h, mu,p0)
% input: h,mu is 2-D data of size p*p.
% p0>=p is the dimension of b
[p1,~]=size(squeeze(h));
h_big=zeros(p0,p0);
p01=fix(p0/2)*2;
p11=fix(p1/2)*2;


ind1=p01/2-p11/2+1; ind2=p01/2+p11/2+1;
h_big(ind1:ind2, ind1:ind2)=squeeze(h).*squeeze(mu);

  gF=fftshift(fft2(ifftshift(h_big)));
end
  