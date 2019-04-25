function BandgapLN(N_AlSb,N_InAs)

% Written by Yajun Wei c2004, updated for HPT by Nate Coirier c2019
% [- (AlSb)m - InSb - (InAs)n - AlAs - ]
% Liquid Nitrogen temperature: 77K
% m layers of AlSb, one layer of InSb, n layers of AlSb, and one layer of AlAs
% m=N_AlSb; n=N_InAs
% One can use parameter R to define the substrate type by matching the ratio, but with more flexibility;
% This model can be used to calculate different interfaces;


global a_AlSb; global a_InSb; global a_InAs; global a_AlAs;
global D_AlSb; global D_InSb; global D_InAs; global D_AlAs;
global a_para; 
global Delta_a_AlSb; global Delta_a_AlAs; global Delta_a_InSb; global Delta_a_InAs; 
global Delta_c_AlSb; global Delta_c_AlAs; global Delta_c_InSb; global Delta_c_InAs; 
global Epa_AlSb; global Epc_AlSb; global Esa_AlSb; global Essa_AlSb; global Esc_AlSb; global Essc_AlSb; global Esasc_AlSb; global Exaxc_AlSb; global Esaxc_AlSb; global Exasc_AlSb; global Exayc_AlSb; global Exassc_AlSb; global Essaxc_AlSb;
global Epa_InSb; global Epc_InSb; global Esa_InSb; global Essa_InSb; global Esc_InSb; global Essc_InSb; global Esasc_InSb; global Exaxc_InSb; global Esaxc_InSb; global Exasc_InSb; global Exayc_InSb; global Exassc_InSb; global Essaxc_InSb;
global Epa_AlAs; global Epc_AlAs; global Esa_AlAs; global Essa_AlAs; global Esc_AlAs; global Essc_AlAs; global Esasc_AlAs; global Exaxc_AlAs; global Esaxc_AlAs; global Exasc_AlAs; global Exayc_AlAs; global Exassc_AlAs; global Essaxc_AlAs;
global Epa_InAs; global Epc_InAs; global Esa_InAs; global Essa_InAs; global Esc_InAs; global Essc_InAs; global Esasc_InAs; global Exaxc_InAs; global Esaxc_InAs; global Exasc_InAs; global Exayc_InAs; global Exassc_InAs; global Essaxc_InAs;
global composition;
global R_InAs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_InAs = 1.0001; % R is the ratio of the lateral superlattice constant to that of the substrate;

Esa_AlSb=-4.55720 - 6.3411; % -6.1714 - 6.5823;
Esc_AlSb=-4.1188 - 6.3411;% -2.0716 - 6.5823;
Epa_AlSb=0.01635 - 6.3411;% 0.9807 - 6.5823;
Epc_AlSb=4.87411 - 6.3411;% 3.0136 - 6.5823;
Essa_AlSb=9.84286 - 6.3411;% 6.7607 - 6.5823;
Essc_AlSb=7.43245 - 6.3411;% 6.1543 - 6.5823;
Esasc_AlSb=-6.63365;% -5.6448;
Esaxc_AlSb=4.58724;% 4.9121;
Exasc_AlSb=8.53398;% 4.2137;
Essaxc_AlSb=7.38446;% 4.3662;
Exassc_AlSb=6.29608;% 3.076;
Exaxc_AlSb=1.10706;% 1.7199;
Exayc_AlSb=4.8996;% 3.6648;
Delta_a_AlSb=0.70373;% 0.973;
Delta_c_AlSb=0.03062;% 0.024;

Esa_InSb=-9.3378 - 5.9071;
Esc_InSb=-3.3248 - 5.9071;
Epa_InSb=0.39352 - 5.9071;
Epc_InSb=2.0791 - 5.9071;
Essa_InSb=6.6378 - 5.9071;
Essc_InSb=5.3807 - 5.9071;
Esasc_InSb=-5.8320;
Esaxc_InSb=4.1129;
Exasc_InSb=7.5769;
Essaxc_InSb=3.4448;
Exassc_InSb=5.8873;
Exaxc_InSb=1.2596;
Exayc_InSb=4.0026;
Delta_a_InSb=0.973;
Delta_c_InSb=0.393;

Esa_InAs = -9.3361 - 6.4851;
Esc_InAs = -3.8093 - 6.4851;
Epa_InAs = 1.159 - 6.4851;
Epc_InAs = 2.8433 - 6.4851;
Essa_InAs=6.8611 - 6.4851;
Essc_InAs=6.0183 - 6.4851;
Esasc_InAs=-6.4183;
Esaxc_InAs=5.0384; 
Exasc_InAs=7.1688;
Essaxc_InAs=3.1749;
Exassc_InAs=5.2785;
Exaxc_InAs=1.9542;
Exayc_InAs=5.0527;
Delta_a_InAs=0.42;
Delta_c_InAs=0.29;

Esa_AlAs=-3.21537- 7.2071; %-7.6201 - 7.2071;
Esc_AlAs=-9.52462- 7.2071; %-1.1823 - 7.2071;
Epa_AlAs=-0.09711- 7.2071; %0.8905 - 7.2071;
Epc_AlAs=4.97139- 7.2071; %3.4939 - 7.2071;
Essa_AlAs=12.05550- 7.2071; %7.3905 - 7.2071;
Essc_AlAs=3.99445- 7.2071; %6.6339 - 7.2071;
Esasc_AlAs=-8.84261; %-6.6642;
Esaxc_AlAs=2.42476; %5.1106;
Exasc_AlAs=13.20317; %6.3;
Essaxc_AlAs=5.83246; %4.58;
Exassc_AlAs=4.60075; %7.3;
Exaxc_AlAs=-0.01434; %1.879;
Exayc_AlAs=4.25494; %3.86;
Delta_a_AlAs=0.29145;
Delta_c_AlAs=0.03152;

Eg_AlSb = 1.68; Eg_InSb = 0.17; Eg_InAs = 0.354; Eg_AlAs = 2.2; % 77K

N = 20 * ( N_AlSb + N_InAs + 2);

% calculating the composition matrix based on the input
% Composition: Al In As Sb
composition = zeros( N/10, 4);
for m = 1:N_AlSb
    composition(2*m-1,1) = 1;
    composition(2*m-1,4) = 1;
end
m = N_AlSb + 1;
composition(2*m,2)=1;
composition(2*m,4)=1;
for m = N_AlSb+2:N_AlSb+N_InAs+1
    composition(2*m-1,2)=1;
    composition(2*m-1,4)= xsegregate(m - (N_AlSb+2));
    composition(2*m-1,3)= 1-composition(2*m-1,4);
end
m = N_AlSb+N_InAs+2;
composition(2*m,1) = 1;
composition(2*m,4)= xsegregate(m - (N_AlSb+2));
composition(2*m,3)= 1-composition(2*m,4);

a_para = 6.089; % SL lattice constant along epi-layer plains

a_AlSb = 6.135; a_InSb = 6.473656; a_InAs = 6.05372*R_InAs;  a_AlAs = 5.66; % Lattice Constants @ 77K

c11_InAs_77K = 8.34; % *1e11 dyn/cm2
c12_InAs_77K = 4.54;
D_InAs = 2*c12_InAs_77K/c11_InAs_77K;

c11_AlSb_77K =8.939;
c12_AlSb_77K =4.427;
D_AlSb = 2*c12_AlSb_77K/c11_AlSb_77K;

c11_AlAs_77K =11.88 + 0.14;
c12_AlAs_77K =5.38 + 0.32;
D_AlAs = 2*c12_AlAs_77K/c11_AlAs_77K;

c11_InSb_77K =6.8;
c12_InSb_77K =3.75;
D_InSb = 2*c12_InSb_77K/c11_InSb_77K;

H0 = Hamiltonian_SL_0(N,N_AlSb,N_InAs);

% find the SL period: az
az = 0;
for m = 1:N_AlSb
   [lattice, tao] = tao_Hac_AlSb(m);
   az = az + lattice(1);
end
m = N_AlSb + 1;
[lattice, tao] = tao_Hac_IF1(m);
az = az + lattice(1); 
for m = N_AlSb+2:N_AlSb+N_InAs+1
   [lattice, tao] = tao_Hac_InAs(m);
   az = az + lattice(1); 
end
m = N_AlSb + N_InAs + 2;
[lattice, tao] = tao_Hac_IF2(m);
az = az + lattice(1);
for m = 1:N_AlSb
   [lattice, tao] = tao_Hca_AlSb(m);
   az = az + lattice(1); 
end
m = N_AlSb + 1;
[lattice, tao] = tao_Hca_IF1(m);
az = az + lattice(1);
for m = N_AlSb+2:N_AlSb+N_InAs+1
    [lattice, tao] = tao_Hca_InAs(m);
    az = az + lattice(1); 
end
m = N_AlSb + N_InAs + 2;
[lattice, tao] = tao_Hca_IF2(m);
az = az + lattice(1); 
az = az/4;

format compact;
a_perp_avr = 2*az / (N_AlSb + N_InAs + 2);
mismatch_LN = (a_perp_avr - a_para)/a_para * 1e6;
[az_RT, mismatch_RT] = Mismatch_RT(N_AlSb,N_InAs);

uz = pi/az;
ux = pi/a_para;
L1 = sqrt(4*ux*ux+uz*uz);
L2 = L1 + sqrt(2)*ux;
L3 = L2 + sqrt(2)*ux;
L4 = L3 + 2*ux;
L5 = L4 + uz;

division = 100;

choice = 'bandgap';
% Specify choice = 'k-space' when you want to calculate the band structures

switch choice

case 'bandgap'
    
    az_RT
    mismatch_RT
    H = Hamiltonian_SL_fg(N,N_AlSb,N_InAs,[0, 0, 0]) .* H0;
    E0 = eig(H);
    E0 = sort(real(E0))';
    Ev_0 = E0(8*(N_AlSb+N_InAs+2))
    Ec_0 = E0(8*(N_AlSb+N_InAs+2)+1);
    Eg = Ec_0-Ev_0
    Wavelength = 1.239/Eg

case 'k-space'
    
    u1 = 1: -1/division: 0;
    Number_of_points_1 = length(u1);
    Energy_1 = zeros(Number_of_points_1,N);
    index = 1;
    for u = 1: -1/division: 0
       E = eig(Hamiltonian_SL_fg(N,N_AlSb,N_InAs,[2*ux*u, 0, uz*u]) .* H0);
       Energy_1(index,:)=E';
       index = index + 1;
    end
    Energy_1 = real(Energy_1);
    Energy_1 = sort(Energy_1')';
    disp('L1 done');
    u2 = 0: 1/division: 1;
    Number_of_points_2 = length(u2);
    Energy_2 = zeros(Number_of_points_2,N);
    index = 1;
    for u = 0: 1/division: 1
       E = eig(Hamiltonian_SL_fg(N,N_AlSb,N_InAs,[ux*u, -ux*u, 0]) .* H0);
       Energy_2(index,:)=E';
       index = index + 1;
    end
    Energy_2 = real(Energy_2);
    Energy_2 = sort(Energy_2')';
    disp('L2 done');
    u3 = 1: 1/division: 2;
    Number_of_points_3 = length(u3);
    Energy_3 = zeros(Number_of_points_3,N);
    index = 1;
    for u = 1: 1/division: 2;
       E = eig(Hamiltonian_SL_fg(N,N_AlSb,N_InAs,[ux*u, ux*(u-2), 0]) .* H0);
       Energy_3(index,:)=E';
       index = index + 1;
    end
    Energy_3 = real(Energy_3);
    Energy_3 = sort(Energy_3')';
    disp('L3 done');
    u4 = 1: -1/division: 0;
    Number_of_points_4 = length(u4);
    Energy_4 = zeros(Number_of_points_4,N);
    index = 1;
    for u = 1: -1/division: 0;
       E = eig(Hamiltonian_SL_fg(N,N_AlSb,N_InAs,[2*ux*u, 0, 0]) .* H0);
       Energy_4(index,:)=E';
       index = index + 1;
    end
    Energy_4 = real(Energy_4);
    Energy_4 = sort(Energy_4')';
    disp('L4 done');
    u5 = 0: 1/division: 1;
    Number_of_points_5 = length(u5);
    Energy_5 = zeros(Number_of_points_5,N);
    index = 1;
    for u = 0: 1/division: 1;
       E = eig(Hamiltonian_SL_fg(N,N_AlSb,N_InAs,[0, 0, uz*u]) .* H0);
       Energy_5(index,:)=E';
       index = index + 1;
    end
    Energy_5 = real(Energy_5);
    Energy_5 = sort(Energy_5')';
    disp('L5 done');

    figure;
    k = [0:L1/division:L1, L1:(L2-L1)/division:L2, L2:(L3-L2)/division:L3, L3:(L4-L3)/division:L4, L4:(L5-L4)/division:L5];
    Energy = [Energy_1; Energy_2; Energy_3; Energy_4; Energy_5];
    plot(k, Energy,'k');
    hold on;
    x1=[L1,L1]; y=[-20,20];
    x2=[L2,L2]; x3=[L3,L3]; x4=[L4,L4];
    plot(x1,y,'g');
    plot(x2,y,'g');
    plot(x3,y,'g');
    plot(x4,y,'g');
    hold off;
    axis([0,L5,-10,-4]);

otherwise
    disp('Wrong Choice! Please check the input value of Choice...');
end