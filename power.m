clear; clc; close all;

np = 10001;

delta = 0.012;
C_D0_c = 0.0093961;
C_D0_t = 2.3856961;
A_R = 19.66;
S = 45.8;
mass = 5500/2.2;
W = mass*9.81;
rho_c = 1.225*0.7422;
rho_t = 1.225;
U_stall = 30.5; 
U_to = U_stall * 1.2;   %must operate above U_to
power_max_t = 236000;
power_max_c = 201000;
D = 2.8;
n = 16;
mu_rolling = 0.03;
C_L_t = 0.75;

U = linspace(0.1, 50, np);
C_L = W ./ (0.5*rho_c*U.^2*S);
C_Di = (1+delta)*C_L.^2 / (pi*A_R);
D_i_c = C_Di * 0.5*rho_c.*U.^2;
D_0_c = C_D0_c*0.5*rho_c.*U.^2;
Drag_c = D_0_c + D_i_c;

[power_c, idx] = min(Drag_c);
U_c = U(idx);

J = U/(D*n);
eff_prop = 0.9717*J.^0.5279;
%eff_prop = -0.85.*(J).*(J-2);

Thrust_t = power_max_t * eff_prop ./ U;

static_thrust = power_max_t*0.0149;     %See Mair p103

for thrust_idx = 1:np
    if U(thrust_idx)<5
        Thrust_t(thrust_idx) = static_thrust;
    end
end

Lift_t = C_L_t*0.5*rho_t.*U.^2*S;
D_0_t = C_D0_t*0.5*rho_c.*U.^2;
D_i_t = ((1+delta)*C_L_t^2 / (pi*A_R))*0.5*rho_c.*U.^2;
Drag_t = D_0_t + D_i_t;

dist(1) = 0;
U_t(1) = 0;
i = 1;
dt = 0.0001;

while U_t(end)<U_to
    [~, j] = min(abs(U-U_t(i)));
    a(i) = (Thrust_t(j) - Drag_t(j) - mu_rolling*(W-Lift_t(j)))/mass;
    U_t(i+1) = U_t(i) + a(i)*dt;
    dist(i+1) = dist(i) + U_t(i)*dt + 0.5*a(i)*dt^2;
    i=i+1;
end

a(end+1)=a(end);
disp([dist(end) "takeoff distance (m)"])
plot(dist, U_t)

%plot(U, Thrust_t)
%hold on;
%plot(U, Drag_t)
%legend(["Thrust", "Drag"])






