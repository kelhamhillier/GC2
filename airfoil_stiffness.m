clear; clc; close all;

np = 10001;

n = 6;

c_r = 1.85;
c_t = 0.65 * c_r
c_av = (c_r+c_t)/2
s = 15;
S = (c_r + c_t)*s
A_R = (4*s^2)/S
W = 5500/2.2*9.81;
U = 42.6;
C_L = 1.41;

rho = 1.225*0.7422;
rho_fuel = 800;
rho_mat = 8000;

thick_sec = 0.15;
thick_mat = 0.015*thick_sec;
m = 0.02;
p = 0.4;

y = linspace(0,s,np);
c = c_r - (c_r - c_t)*(y/s);

x = linspace(0,1,np);
t = 10*thick_sec*(0.2969*x.^0.5 - 0.126*x -0.3516*x.^2 + 0.2843*x.^3 - 0.1015*x.^4);
for i = 1:np
    if x(i)<=p
        y_c(i) = (m/p^2)*(2*p*x(i) - x(i)^2);
    else
        y_c(i) = (m/(1-p)^2)*(1-2*p+2*p*x(i)-x(i)^2);
    end
end

A = sum(t)/(np);
x_g = sum(x.*t)/((np)*A);
y_g = sum(y_c.*t)/((np)*A);

I_xx = sum(((t/2)+(y_c-y_g)).^2 * thick_mat)/(np) + sum(((t/2)-(y_c-y_g)).^2 * thick_mat)/(np);

for j = 1:np
    A_c(j) = A*c(j)^2;  %Area of section
    mass(j) = A_c(j)*rho_fuel + thick_mat*c(j)*2*rho_mat;  %Mass per unit length
end

M(np) = 0;
Shear(np) = 0;

Gamma_0 = W / (rho * U * (pi/2) * s);
lambda = (pi*Gamma_0/4)/(1 - (1 - c_t/c_r)/2);

L = n*rho*U/2 * (Gamma_0 * (1 - y.^2/s^2).^0.5 + lambda*(1-(1-c_t/c_r)*y/s));    %Loading per unit length
dL = L * (s/np);
dW = n*mass*9.81*(s/np);

U_de = 66 * 0.3048;
mu = (2*W/S) / (rho * 9.81 * c_av * C_L);
K = (0.88*mu) / (5.3 + mu);
U_gust = K * U_de;
n_gust = (rho*U_gust*U*C_L)/(2*W/S);

for k = 1:np-1
    Shear(np-k) = Shear(np-k+1) + dW(np-k+1) - dL(np-k+1);
    M(np-k) = M(np-k+1) + (Shear(np-k+1))*(s/np);
end

disp(["Moment at root", M(1)])

max_stress = M(1) * max(abs(t/2 + y_g - y_c))*c(1) / (I_xx*c(1)^4);

disp(["Maximum stress " max_stress/10^6 "MPa"])

disp([sum(A_c)*s/np  "volume per wing (m^3)"]);

disp([sum(A_c)*s/np * 2 * rho_fuel "kg of wing fuel"])

disp([sum(mass)*s/np * 2 "total mass of wings (kg)"])

plot(y, L./(0.5*rho*U^2*c));
xlabel("y")
ylabel("c_l")
%hold on;
%plot(y, n*rho*U * (Gamma_0 * (1 - y.^2/s^2).^0.5));
%plot(y, n*rho*U * lambda*(1-(1-c_t/c_r)*y/s));
