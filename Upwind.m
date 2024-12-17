function [f] = Upwind(fj_1, fj, fj1)
%Upwind function is used to build 3th order upwind (positive vector)
f = fj+0.5*(fj1-fj_1);
end

