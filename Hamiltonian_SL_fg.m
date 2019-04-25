function matrix = Hamiltonian_SL_fg(N,N_AlSb,N_InAs,kvector)

% Written by Yajun Wei c2004

% This function is to generate the superlattice Hamiltonian
% The model takes into account: nearest neighbor interaction, strain effects.

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


matrix = ones(N);

for m = 1:N_AlSb
   [lattice, tao] = tao_Hac_AlSb(m);
   matrix( 20*m-19:20*m-10, 20*m-9:20*m ) = Hac_fg(kvector,tao);
end
m = N_AlSb + 1;
[lattice, tao] = tao_Hac_IF1(m);
matrix( 20*m-19:20*m-10, 20*m-9:20*m ) = Hac_fg(kvector,tao); 
for m = N_AlSb+2:N_AlSb+N_InAs+1
   [lattice, tao] = tao_Hac_InAs(m);
   matrix( 20*m-19:20*m-10, 20*m-9:20*m ) = Hac_fg(kvector,tao); 
end
m = N_AlSb + N_InAs + 2;
[lattice, tao] = tao_Hac_IF2(m);
matrix( 20*m-19:20*m-10, 20*m-9:20*m ) = Hac_fg(kvector,tao);

for m = 1:N_AlSb
   [lattice, tao] = tao_Hca_AlSb(m);
   matrix( 20*m-9:20*m, 20*m+1:20*m+10 ) = Hca_fg(kvector,tao); 
end
m = N_AlSb + 1;
[lattice, tao] = tao_Hca_IF1(m);
matrix( 20*m-9:20*m, 20*m+1:20*m+10 ) = Hca_fg(kvector,tao);
for m = N_AlSb+2:N_AlSb+N_InAs+1
    [lattice, tao] = tao_Hca_InAs(m);
    matrix( 20*m-9:20*m, 20*m+1:20*m+10 ) = Hca_fg(kvector,tao); 
end

m = N_AlSb + N_InAs + 2;
[lattice, tao] = tao_Hca_IF2(m);
matrix( 1:10, 20*m-9:20*m ) = Hca_fg(kvector,tao)'; 

matrix = matrix .* matrix';

