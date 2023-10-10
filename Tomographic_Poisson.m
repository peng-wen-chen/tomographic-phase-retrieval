% use CG to compute f from h_k
% In this code, we conduct the reconstruction
% of tomographic phase retrieval under  one random mask 
% The object can be a phantom or a RPP. 
clear
input.RPP=0; % =1, if it is RPP; =0 if it is phantom
input.num=1; % # of random masks
input.NSR=0.00; %input.NSR=0.6;
n=36; % n should be a square.
p=2*n-1; % p (the size of h is p*p) should be at least 2n-1
p0=p; %(the size of g is p0*p0)
m=n;
input.n1=n^(3/2);
input.s=fix(n/2)+1; % the location, where mask plane is placed 
%%
input.n=n;
input.p=p;
input.p0=p0;

input.m=m;
input.CGiter=10; % # of iterations in  inexact CG
input.Nuliter=10;  % # of power iterations in  null vector method
input.APiter=10;  % # of iterations in AP

%% set up the index range for support
p1=fix(p/2)*2; % p1 is the greatest even integer, not exceeding p  
input.p1=p1;
% 
% the central index interval[index1, index2] is  the support interval.
index1=p1/2-fix(n/2)+1;
input.ind1=index1; 
index2=p1/2-fix(n/2)+n;
input.ind2=index2; 
input.rho=[-p1/2:p1/2]; % Z_p
%%  diffraction patterns
input = Tomo_initial(input); 
b1=input.b1;
b2=input.b2;
b3=input.b3;
% generate indicators 
 weight=weight_gen(input);
 % initialization 
fs=input.fs;
   f0out=threeD_to_twoD(fs(index1: index2,index1: index2,index1: index2 ),input);
figure; 
imshow(f0out);  % show the real part of the unknown object

%% Simulation starts with random initialization fs
rng('default')
rng(1)
f1=randn(p,p,p)+i*randn(p,p,p); 
fs=zeros(p,p,p);

fs(index1:index2, index1:index2,index1:index2)=f1(index1: index2,index1: index2,index1: index2); % zero padding image.
fs=fs/norm(fs(:))*norm(input.fs(:));

 norm_f00=norm(input.fs(:));
 fs=fs/norm(fs(:))* norm_f00;
 null_ini=fs;
 
 %% Null vector call
  for kk=1:1
      [null,cor] = null_vector(weight,null_ini, input);
           null_ini=null;
  end
   null2D=threeD_to_twoD(null(index1: index2,index1: index2,index1: index2),input);
figure(12); % show the moduli of null vector
imshow(abs(null2D))
   %%   AP one axis function call
   b.data1=input.b1;
   b.data2=input.b2;
   b.data3=input.b3;
   fs_ini=null;
 for kkk=1:1  
[fs, err]=AP_method(b,fs_ini,input); fs_ini=fs;
f0out=threeD_to_twoD(fs(index1:index2,index1: index2,index1: index2),input);
figure(1); 
imshow(abs(f0out))
 end
 
 
 