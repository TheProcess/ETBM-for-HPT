function [az, mismatch] = Mismatch_RT(N_AlSb,N_InAs)

% Written by Yajun Wei c2004
% [- (AsIn)m - As - Gax1In(1-x1) - (SbGa)n - Sb - Gax2In(1-x2) - ] N
% Liquid Nitrogen temperature: 77K
% m layers of AsIn, one layer of AsGaxIn(1-x), n layers of SbGa, and one layer of SbGaxIn(1-x)
% m=N_AlSb; n=N_InAs
% R is the ratio of the lateral superlattice constant to that of the substrate;
% One can use parameter R to define the substrate type by matching the ratio, but with more flexibility;
% This model can be used to calculate different interfaces;

global a_AlSb; global a_InSb; global a_InAs; global a_AlAs;
global D_AlSb; global D_InSb; global D_InAs; global D_AlAs;
global a_para; 
global Epa_AlSb; global Epc_AlSb; global Esa_AlSb; global Essa_AlSb; global Esc_AlSb; global Essc_AlSb; global Esasc_AlSb; global Exaxc_AlSb; global Esaxc_AlSb; global Exasc_AlSb; global Exayc_AlSb; global Exassc_AlSb; global Essaxc_AlSb;
global Epa_InSb; global Epc_InSb; global Esa_InSb; global Essa_InSb; global Esc_InSb; global Essc_InSb; global Esasc_InSb; global Exaxc_InSb; global Esaxc_InSb; global Exasc_InSb; global Exayc_InSb; global Exassc_InSb; global Essaxc_InSb;
global Epa_AlAs; global Epc_AlAs; global Esa_AlAs; global Essa_AlAs; global Esc_AlAs; global Essc_AlAs; global Esasc_AlAs; global Exaxc_AlAs; global Esaxc_AlAs; global Exasc_AlAs; global Exayc_AlAs; global Exassc_AlAs; global Essaxc_AlAs;
global Epa_InAs; global Epc_InAs; global Esa_InAs; global Essa_InAs; global Esc_InAs; global Essc_InAs; global Esasc_InAs; global Exaxc_InAs; global Esaxc_InAs; global Exasc_InAs; global Exayc_InAs; global Exassc_InAs; global Essaxc_InAs;
global composition;
global Delta_a_AlSb; global Delta_a_AlAs; global Delta_a_InSb; global Delta_a_InAs; 
global Delta_c_AlSb; global Delta_c_AlAs; global Delta_c_InSb; global Delta_c_InAs; 
global R_InAs;

a_para = 6.09593; % SL lattice constant along epi-layer plains

a_AlSb = 6.0583; a_InSb = 6.47937; a_InAs = 6.09593*R_InAs;  a_AlAs = 5.65325; % Lattice Constants @ room temperature

D_AlSb = 1.086805139; D_InSb = 1.093117409; D_InAs = 0.910799185;  D_AlAs = 0.903610411; % Strain D001 Constant @ room temperature

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

a_perp_avr = 2*az / (N_AlSb + N_InAs + 2);
mismatch = (a_perp_avr - a_para)/a_para * 1e6;

%restore the values for low temperature calculations
a_AlSb = 6.05372; a_InSb = 6.473656; a_InAs = 6.089*R_InAs;  a_AlAs = 5.648206; % Lattice Constants @ 77K

c11_InAs_77K = 9.05; % *1e11dyn/cm2
c12_InAs_77K = 4.115;
D_InAs = 2*c12_InAs_77K/c11_InAs_77K;

c11_AlSb_77K =8.3;
c12_AlSb_77K =4.95;
D_AlSb = 2*c12_AlSb_77K/c11_AlSb_77K;

c11_AlAs_77K =12.059;
c12_AlAs_77K =5.41072;
D_AlAs = 2*c12_AlAs_77K/c11_AlAs_77K;

c11_InSb_77K =6.8;
c12_InSb_77K =3.75;
D_InSb = 2*c12_InSb_77K/c11_InSb_77K;
