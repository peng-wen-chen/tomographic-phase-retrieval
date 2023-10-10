function fout=threeD_to_twoD(f1,input)
% f1 is 3-D object with size n*n*n
% fout is 2-D object with size n^(3/2) * n^(3/2)
n=input.n;
for kk=1:n
i2=rem(kk-1,n^0.5)+1;
j2=ceil(kk/(n^0.5));
fout((i2-1)*n+1:i2*n,(j2-1)*n+1:j2*n)=f1(kk,1:n,1:n);
end
end