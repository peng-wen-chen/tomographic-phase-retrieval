function fs = P_omega(fs_ini, input)
% Apply P_{Z_n}
p=input.p;
fs=zeros(p,p,p);
ind1=input.ind1; 
ind2=input.ind2; 
fs(ind1:ind2, ind1:ind2, ind1:ind2)=fs_ini(ind1:ind2, ind1:ind2, ind1:ind2);

if input.RPP==0
    
fs=max(real(fs),0);
else
end


end

