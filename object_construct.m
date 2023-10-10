function [fp,fcube] = object_construct(input)
% generate phantom or RPP (according to input.RPP) as 
% Output: fp is 2D object, fcube is 3D object.
n1=input.n1;
RPPmask=exp(i*2*pi*rand(n1,n1)*input.RPP); %  RPP or phantom
fp=phantom(n1,n1).*RPPmask;
fcube=twoD_to_threeD(fp,input);
end

