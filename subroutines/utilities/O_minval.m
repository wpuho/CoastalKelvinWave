function [mx]=O_minval(a)

%eda:show the minimum value of the arrary
[m n l k]=size(a);

if l==1 %2D

mx=min(min(a));

elseif l~=1 & k==1 % 3D
mx=min(min(min(a)));

elseif l~=1 & k~=1 %4D
mx=min(min(min(min(a))));

end

return;
