function    weight=weight_gen(input)
   % construct the indicators for null vector methods
   % compute the norm of diffraction patterns
   m=input.m;
   % x-lines
for ii=1:m
    tt=input.b1(ii,:,:);
    ttt=sort(tt(:));
    mean_b=ttt(fix(end/2));
ind1(ii,:,:)=(tt>mean_b);
% compute the norm of diffraction patterns
ind1_norm(ii)=norm(squeeze(input.b1(ii,:,:)),'fro');
end
  % y-lines
for ii=1:m
    tt=input.b2(:,ii,:);
    ttt=sort(tt(:));
    mean_b=ttt(fix(end/2));
ind2(:,ii,:)=(tt>mean_b);
ind2_norm(ii)=norm(squeeze(input.b2(:,ii,:)),'fro');
end
% z-lines
for ii=1:m
    tt=input.b3(:,:,ii);
    ttt=sort(tt(:));
    mean_b=ttt(fix(end/2));
 ind3(:,:,ii)=(tt>mean_b);
ind3_norm(ii)=norm(squeeze(input.b3(:,:,ii)),'fro');
end
   
   weight.data1=ind1;
   weight.data2=ind2;
   weight.data3=ind3;
   weight.norm1=ind1_norm;
   weight.norm2=ind2_norm;
   weight.norm3=ind3_norm;
end