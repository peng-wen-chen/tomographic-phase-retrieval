function LHS=CG_LHS(fs,T_hat, p,n)
% the code computes LHS of CG normal equation.
p1=floor(p/2)*2;
% F_fs is the Fourier transform F{f}
F_fs=fftshift(fftn(ifftshift(fs)))/p^(3/2);

%
FF=(p^0.5)*T_hat.*F_fs;
FF_ifft=fftshift(ifftn(ifftshift(FF)))*p^(3/2);
  % Apply the support projection P_{Z_n} 
 LHS=zeros(p,p,p);
  ind1=p1/2+1-fix(n/2);
ind2=p1/2+n-fix(n/2);
 LHS(ind1:ind2,ind1:ind2,ind1:ind2)=FF_ifft(ind1:ind2,ind1:ind2,ind1:ind2);

end
