function f0=twoD_to_threeD(fp,input)
% Construct 3-D object   f0 from 2-D object fp
% by stacking the slices from partitioning fp
n=input.n;

f0=zeros(n,n,n);
for kk=1:n
i2=rem(kk-1,n^0.5)+1;
j2=ceil(kk/(n^0.5));

f0(kk,1:n,1:n)=fp((i2-1)*n+1:i2*n,(j2-1)*n+1:j2*n);
end

end