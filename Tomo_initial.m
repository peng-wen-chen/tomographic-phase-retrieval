function input = Tomo_initial(input)

if input.NSR>1e-10 % if less than 10^(-10), it is regarded noise-free
ss=1/input.NSR; 
else
    input.noise1=0;
ss=0;    
end

n=input.n;
p=input.p;
p0=input.p0;

m=input.m;

%% generate m dirctionvectors, 
%  beta, two of whose components randomly in [-1,1]
rand('state',0)
 beta1(1:input.num:m,1:2)=rand(m/input.num,2).*(2*(rand(m/input.num,2)>0.5)-1);
 beta2(1:input.num:m,1:2)=rand(m/input.num,2).*(2*(rand(m/input.num,2)>0.5)-1);
 beta3(1:input.num:m,1:2)=rand(m/input.num,2).*(2*(rand(m/input.num,2)>0.5)-1);
 for ii=2:input.num
beta1(ii:input.num:m,1:2)=beta1(1:input.num:m,1:2);
beta2(ii:input.num:m,1:2)=beta2(1:input.num:m,1:2);
beta3(ii:input.num:m,1:2)=beta3(1:input.num:m,1:2);
end
input.beta1=beta1;
input.beta2=beta2;
input.beta3=beta3;

% generate random masks
for kk=1:m
rand('state',rem(kk,input.num));    
%mu=exp(i*2*pi*rand(p,p)); 
input.mu(kk,1:p,1:p)=exp(i*2*pi*rand(p,p));
end
% generate object (phantom or RPP)
[fp,f0]=object_construct(input);


% %%
%    if rem(p,2)==1
%           p1=p-1;
%     else
%         p1=p;
%     end
% 
%% generate f0
f=zeros(p,p,p);
index1=input.ind1; 
index2=input.ind2; 
f(index1:index2, index1:index2,index1:index2)=f0; % zero padding image.
input.fs=f;
input.fp=fp;
%%
 % construct kernel T_hat
input.T=T_construct(input);
%% construct diffraction patterns b1,b2,b3
% R-map
[fb1,fb2,fb3]=Rmap(input.fs, input);
gF1=zeros(m,p0,p0);
gF2=zeros(p0,m,p0);
gF3=zeros(p0,p0,m);

% Q-map
for ii=1:m
gF1(ii,:,:)=Qmap(fb1(ii,:,:), squeeze(input.mu(ii,:,:)),p0);
gF2(:,ii,:)=Qmap(fb2(:,ii,:), squeeze(input.mu(ii,:,:)),p0);
gF3(:,:,ii)=Qmap(fb3(:,:,ii), squeeze(input.mu(ii,:,:)),p0);
 end

%%  Add noise to  noise-free patterns b1,b2,b3
    b1=abs(gF1); %input.b1=b1;
        b2=abs(gF2);% input.b2=b2;
    b3=abs(gF3);% input.b2=b2;

        
        input.noise1=(norm([b1(:).^0.5;b2(:).^0.5;b3(:).^0.5])...
            /norm([b1(:); b2(:); b3(:)]))^2*ss *1e-6;
  
    input.b1=zeros(size(gF1));
    input.b2=zeros(size(gF2));
    input.b3=zeros(size(gF3));
    
    
    if input.noise1>0
    input.b1=(imnoise((input.noise1*b1).^2,'poisson')).^0.5/(input.noise1);
    else
        % noise-free case
        input.b1=b1;
    end
 
    
    if input.noise1>0
    input.b2=(imnoise((input.noise1*b2).^2,'poisson')).^0.5/(input.noise1);
    else
        input.b2=b2;
    end
    
    if input.noise1>0
            input.b3=(imnoise( (input.noise1*b3).^2,'poisson')).^0.5/(input.noise1);
    else
        input.b3=b3;
    end


    
% disp('Poisson noise level is ')
% norm(input.b1(:).^2-b1(:).^2)^0.5/norm(b1(:))
% norm(input.b2(:).^2-b2(:).^2)^0.5/norm(b2(:))
% norm(input.b3(:).^2-b3(:).^2)^0.5/norm(b3(:))

%     
% input.u1=gF1./b1;
% input.u2=gF2./b2;
% input.u3=gF3./b3;

end 