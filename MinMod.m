function [f] = MinMod(fj_1,fj,fj1)
%MinMod function is used to create Min-mod limitor for the space discretilaztion
r=zeros(1,3);
phi=zeros(1,3);
for i=1:3
    if fj1(i)==fj(i)  && fj_1(i)==fj(i)
        r(i)=1;
    % elseif  fj1(i)==fj(i)  && fj_1(i)~=fj(i)
    %     r(i)=-1;
    % elseif  fj1(i)~=fj(i)   && fj_1(i)==fj(i)
    %     r(i)=-1 ;
    else
        r(i)=(fj(i)-fj_1(i))/(fj1(i)-fj(i));
    end
    if r(i)<0
        phi(i)=0;
    else
        phi(i)=min([1,r(i)]);
    end
    
end
f = fj + 0.5*phi.*(fj1-fj);
end

