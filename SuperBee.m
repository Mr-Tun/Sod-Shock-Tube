function [f] = SuperBee(fj_1,fj,fj1)
%Superbee is used to create superbee limiter for U_(i+0.5)
r=zeros(1,3);
phi=zeros(1,3);
for i=1:3
    if fj1(i)==fj(i)
    r(i)=1;
    else
        r(i)=(fj(i)-fj_1(i))/(fj1(i)-fj(i));
    end
        
    phi(i)=max([0,min(2*r(i),1),min(r(i),2)]);
    

end
f = fj + 0.5*phi.*(fj1-fj);
end

