function [G delta] = induced(theta, gamma, U, s, np)

ng = 1500;
G = zeros(1,ng);
delta = 0;

for i = 1:ng
    integral =0;
    for j = 2:np
        integral = integral+((gamma(j)+gamma(j-1))/2)*(theta(j)-theta(j-1))*sin(i*(theta(j-1)+theta(j))/2);
    end
    G(i) = integral*(2/pi)/(U*s);
    delta = delta + i*(G(i)/G(1))^2;
end
delta = delta-1;
