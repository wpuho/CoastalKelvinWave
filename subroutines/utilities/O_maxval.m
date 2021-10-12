function [mx]=O_maxval(a)

%eda:show the maximum value of the arrary
[m n l k]=size(a);

if l==1 %2D

mx=max(max(a));

elseif l~=1 & k==1 % 3D
mx=max(max(max(a)));

elseif l~=1 & k~=1 %4D
mx=max(max(max(max(a))));

end

return;
