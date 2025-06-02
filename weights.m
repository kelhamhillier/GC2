clear; clc; close all;

m_ft = 3.2808399;   %meters to feet
m_in = 39.3700787;  %meters to inches
kg_lbs = 2.20462262;

A_R = 19.7;
S_w = 45.8 * m_ft^2;   %Wing area
W_fw = 728.5 * kg_lbs;    %Weight of fuel in wing
rho = 0.7422 * 1.225;
U = 36.6;
a = sqrt(1.4 * 287 * 0.9324 * 288.15);
q = 0.5 * rho * U^2 * 0.020885;     %lb/ft^2
lambda = 0.65;      %Taper ratio
lambda_h = 1;
lambda_v = 1;
Lambda = 0;
Lambda_ht = 0;
Lambda_vt = 0;
c = 1.53;       %Average chord of wing
t = c*0.15;       %Average thickness of wing
N_z = 3.8*1.5;     %Ultimate load factor (?)
W_dg = 5500;    %Design gross weight
S_ht = 1.15 * 6.9 * m_ft^2;    %Horizontal tail area
S_vt = 1.15 * 1.5 * m_ft^2;    %Vertical tail area
S_f = 25.5 * m_ft^2;     %Fuselage wetted area
L_t = 6.9 * m_ft;     %Tail length
L = 9.3 * m_ft;       %Fuselage structural length
D = 1 * m_ft;       %Fuselage structural depth
N_l = 0.05;     %Ultimate landing load factor
W_l = 3190;     %Landing design gross weight
L_m = 1.5 * m_in;     %Length of main landing gear
L_n = 1.5 * m_in;     %Nose gear length
W_en = 138.5;   %Engine weight
N_en = 1;   %Number of engines
V_t = 29.8;     %Total fuel volume
V_i = V_t;     %Integral tanks volume
N_t = 2;     %Number of fuel tanks
B_w = 30 * m_ft;    %Wing span
W_uav = 500;       %Uninstalled avionics weight, typically 800-1400lbs
N_p = 2;        %Number of pilots
M = U/a;

W_pilots = 180;

wing_cg = 3.3;
horiz_cg = wing_cg + 6;
vert_cg = horiz_cg;
fuse_cg = 4.3;
main_lg_cg = 6.5-0.9;
nose_lg_cg = 2;
engine_cg = 9.485-0.9;
controls_cg = 0.5;
avionics_cg = 1;
pilots_cg = 1.5;
payload_cg = 3;

W_wing = (0.036* S_w^0.758 * W_fw^0.0035 * (A_R / (cos(Lambda)^2))^0.6 * q^0.006 * lambda^0.04 * ((100*t)/(c*cos(Lambda)))^(-0.3) * (N_z * W_dg)^0.49)/2.2
W_horiztail = (0.016 * (N_z * W_dg)^0.414 * q^0.168 * S_ht^0.896 * ((100*t)/(c*cos(Lambda_ht)))^(-0.12) * (A_R / (cos(Lambda_ht)^2))^0.043 * lambda_h^(-0.02))/2.2
W_verttail = (0.073 * (N_z * W_dg)^0.376 * q^0.122 * S_vt^0.873 * ((100*t)/(c*cos(Lambda_vt)))^(-0.49) * (A_R / (cos(Lambda_vt)^2))^0.357 * lambda_v^0.039)/2.2
W_fuselage = (0.052*S_f^1.086 * (N_z * W_dg)^0.177 * L_t^(-0.051) * (L/D)^(-0.072) * q^0.241)/2.2 
W_maingear = (0.095 * (N_l * W_l)^(0.768) * (L_m / 12)^0.409)/2.2
W_nosegear = (0.125 * (N_l * W_l)^0.566 * (L_n / 12)^0.845)/2.2
W_engine = (2.575 * W_en^0.922 * N_en)/2.2
W_fuelsys = (2.49 * V_t^0.726 * (1/(1 + V_i/V_t))^0.363 * N_t^0.242 * N_en^0.157)/2.2
W_controls = (0.053 * L^1.536 * B_w^0.371 * (N_z * W_dg * 10^(-4))^0.8)/2.2
W_hydraulics = (0.001 * W_dg)/2.2
W_avionics = (W_uav)/2.2
W_elec = (12.57 * (W_avionics+W_fuelsys) ^ 0.51)/2.2

W_empty = W_fuelsys + W_wing + W_horiztail + W_verttail + W_fuselage + W_maingear + W_nosegear + W_engine + W_controls + W_hydraulics + W_avionics + W_elec
W_no_fuel = W_empty + W_pilots + 227
W_fuel = (W_dg - 2*170 - 500) / 2.2 - W_empty
V_t = W_fuel / 800 * m_ft^3

%W_fuel=0;
cg = (W_fuel*(wing_cg-(0.5-0.17)*1.53)+W_wing*wing_cg + W_horiztail*horiz_cg + W_verttail*vert_cg + W_fuselage*fuse_cg + W_maingear*main_lg_cg + W_nosegear*nose_lg_cg + W_engine*engine_cg + W_controls*controls_cg + 0.5*W_avionics*avionics_cg + W_pilots*pilots_cg + 227*payload_cg)/(W_no_fuel+W_fuel)
x_g=cg-wing_cg+(0.25*1.53)
MAC = (cg-wing_cg)/1.53+0.5