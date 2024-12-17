function [f] = VanLeer(fj_1,fj,fj1)
%VanLeer function is used to create Van-Leer limiter for U_(i+0.5)
r=zeros(1,3);
phi=zeros(1,3);
for i=1:3
    if fj1(i)==fj(i)
    r(i)=1;
    else
        r(i)=(fj(i)-fj_1(i))/(fj1(i)-fj(i));
    end
phi(i)=(r(i)+abs(r(i)))/(1+abs(r(i)));
end
f = fj + 0.5*phi.*(fj1-fj);
end

