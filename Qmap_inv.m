function h=Qmap_inv(gF, mu,input)
g=ifftshift(squeeze(gF));
gF=fftshift(ifft2(g));
p0=input.p0;
p1=input.p1;
p01=fix(p0/2)*2;
ind1=p01/2-p1/2+1; ind2=p01/2+p1/2+1;
h=gF(ind1: ind2, ind1:ind2).*conj(mu);
end
  