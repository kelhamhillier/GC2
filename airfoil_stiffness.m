clear; clc; close all;

np = 10001;

n = 6;

c_r = 1.85;
c_t = 0.65 * c_r;
c_av = (c_r+c_t)/2;
s = 15;
S = (c_r + c_t)*s;
A_R = (4*s^2)/S;
W = (5500/2.2)*9.81;
U = 41;
C_L = 1.41;

rho = 1.225*0.7422;
rho_fuel = 800;
rho_mat = 2810;
E = 71.7 * 10^9; 

mass_fuel_total = 489;

thick_sec = 0.15;
thick_mat = 0.011*thick_sec;
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
x_g_fuel = sum(x.*t)/((np)*A);
y_g_fuel = sum(y_c.*t)/((np)*A);

A_metal = thick_mat*2;
x_g = sum(x.*thick_mat*2)/((np)*A_metal);
y_g = sum(y_c.*thick_mat*2)/((np)*A_metal);

I_xx = sum(((t/2)+(y_c-y_g)).^2 * thick_mat)/(np) + sum(((t/2)-(y_c-y_g)).^2 * thick_mat)/(np);

mass_fuel = zeros(1,np);

for j = 1:np
    A_c(j) = A*c(j)^2;  %Area of section
    mass_frame(j) = thick_mat*c(j)*2*rho_mat;  %Mass per unit length
end

for l = 1:np
    if sum(mass_fuel*s/np)<(mass_fuel_total/2)
        mass_fuel(np-l+1) = A_c(np-l+1)*rho_fuel;
    end
end

mass = mass_frame + mass_fuel;

M(np) = 0;
Shear(np) = 0;
Torsion(np) = 0;

Gamma_0 = W / (rho * U * (pi/2) * s);
lambda = (pi*Gamma_0/4)/(1 - (1 - c_t/c_r)/2);

L = n*rho*U/2 * (Gamma_0 * (1 - y.^2/s^2).^0.5 + lambda*(1-(1-c_t/c_r)*y/s));    %Loading per unit length
dL = L * (s/np);
dW = n*mass*9.81*(s/np);

T = mass_fuel*9.81*n*(x_g_fuel-x_g) .*c;
dT = T * (s/np);

%For delta but does not work
theta = linspace(0, pi, np);
y_theta = -s * cos(theta);
gamma = 1/2 * (Gamma_0 * (1 - y_theta.^2/s^2).^0.5 + lambda*(1-(1-c_t/c_r)*abs(y_theta)/s));    %Loading per unit length
[G delta] = induced(theta, gamma, U, s, np); 

U_de = 66 * 0.3048;
mu = (2*W/S) / (rho * 9.81 * c_av * C_L);
K = (0.88*mu) / (5.3 + mu);
U_gust = K * U_de;
n_gust = (rho*U_gust*U*C_L)/(2*W/S);

for k = 1:np-1
    Shear(np-k) = Shear(np-k+1) + dW(np-k+1) - dL(np-k+1);
    M(np-k) = M(np-k+1) + (Shear(np-k+1))*(s/np);
    Torsion(np-k) = Torsion(np-k+1) - dT(np-k+1) + dL(np-k+1)*(x_g - 0.25)*(c(np-k+1)+c(np-k))/2;
end

kappa = M ./ (I_xx * E * c.^4);
psi(1) = 0;
for m = 2:np
    psi(m) = psi(m-1) + kappa(m-1) * (x(m)-x(m-1));
end
dis(1) = 0;
for n = 2:np
    dis(n) = dis(n-1) + psi(n-1) * (x(n)-x(n-1));
end

disp(["Displacement at root", dis(np)*1000, "mm"])

disp(["Moment at root", M(1)])

disp(["Torsion at root", Torsion(1)])

[max_stress loc] = max(((2*abs(Torsion/(2*A.*c.^2*thick_mat))).^2+abs(M .* max(abs(t/2 + y_g - y_c).*c ./ (I_xx.*c.^4))).^2).^0.5);

disp(["Maximum stress" max_stress/10^6 "MPa"])

disp(["at y =" y(loc)])

disp([sum(A_c)*s/np  "volume per wing (m^3)"]);

disp([sum(mass_frame)*s/np * 2 "frame mass of wings (kg)"])

disp([sum(mass)*s/np * 2 "total mass of wings (kg)"])
if n==1
    plot(y, L./(0.5*rho*U^2*c));
    xlabel("y")
    ylabel("c_l")
end
%plot(y, dL);
%hold on;
%plot(y, dW);
