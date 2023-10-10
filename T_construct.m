function [T_hat] = T_construct(input)
% construct T_hat
n=input.n;
p=input.p;
p1=input.p1;
rho=input.rho; 
beta1=input.beta1;
beta2=input.beta2;
beta3=input.beta3;


T1=zeros(p,p,p);
T2=zeros(p,p,p);
T3=zeros(p,p,p);


 for jj=1:p
for kk=1:p 
    % x-lines
 t_ker=sum(exp(-i*2*pi*[0:p1/2]'*[rho(jj),rho(kk)]*beta1'/p),2);%/p^2;
 t_kernel=[t_ker;   conj(t_ker(end:-1:2)) ];
 T1(:,jj,kk)=fftshift(fft(t_kernel,p))/p^0.5;
%T1(:,jj,kk)=fftshift(temp);

% y-lines
 t_ker=sum(exp(-i*2*pi*[0:p1/2]'*[rho(jj),rho(kk)]*beta2'/p),2);
 t_kernel=[t_ker;   conj(t_ker(end:-1:2)) ];
 T2(jj,:,kk)=fftshift(fft(t_kernel,p))/p^0.5;
%T2(jj,:,kk)=fftshift(temp);

% z-lines
 t_ker=sum(exp(-i*2*pi*[0:p1/2]'*[rho(jj),rho(kk)]*beta3'/p),2);
 t_kernel=[t_ker;   conj(t_ker(end:-1:2)) ];
 T3(jj,kk,:)=fftshift(fft(t_kernel,p))/p^0.5;
%T3(jj,kk,:)=fftshift(temp);

end
 end
 T_hat=T1+T2+T3;

end

