function matrix = Hamiltonian_SL_0(N,N_AlSb,N_InAs)

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


matrix = zeros(N);


for m = 1:N_AlSb
   [lattice, E_vec, tao] = parametergen_Hac_AlSb(m);
   matrix( 20*m-19:20*m-10, 20*m-9:20*m ) = Hac(1, 1, lattice(5), lattice(2), lattice(3), lattice(4), E_vec(9), E_vec(10), E_vec(11), E_vec(14), E_vec(15), E_vec(13), E_vec(12)); 
end
m = N_AlSb + 1;
[lattice, E_vec, tao] = parametergen_Hac_IF1(m);
matrix( 20*m-19:20*m-10, 20*m-9:20*m ) = Hac(1, 1, lattice(5), lattice(2), lattice(3), lattice(4), E_vec(9), E_vec(10), E_vec(11), E_vec(14), E_vec(15), E_vec(13), E_vec(12)); 
for m = N_AlSb+2:N_AlSb+N_InAs+1
   [lattice, E_vec, tao] = parametergen_Hac_InAs(m);
   matrix( 20*m-19:20*m-10, 20*m-9:20*m ) = Hac(1, 1, lattice(5), lattice(2), lattice(3), lattice(4), E_vec(9), E_vec(10), E_vec(11), E_vec(14), E_vec(15), E_vec(13), E_vec(12)); 
end
m = N_AlSb + N_InAs + 2;
[lattice, E_vec, tao] = parametergen_Hac_IF2(m);
matrix( 20*m-19:20*m-10, 20*m-9:20*m ) = Hac(1, 1, lattice(5), lattice(2), lattice(3), lattice(4), E_vec(9), E_vec(10), E_vec(11), E_vec(14), E_vec(15), E_vec(13), E_vec(12));

for m = 1:N_AlSb
   [lattice, E_vec, tao] = parametergen_Hca_AlSb(m);
   matrix( 20*m-9:20*m, 20*m+1:20*m+10 ) = Hca(1, 1, lattice(5), lattice(2), lattice(3), lattice(4), E_vec(9), E_vec(10), E_vec(11), E_vec(14), E_vec(15), E_vec(13), E_vec(12)); 
end
m = N_AlSb + 1;
[lattice, E_vec, tao] = parametergen_Hca_IF1(m);
matrix( 20*m-9:20*m, 20*m+1:20*m+10 ) = Hca(1, 1, lattice(5), lattice(2), lattice(3), lattice(4), E_vec(9), E_vec(10), E_vec(11), E_vec(14), E_vec(15), E_vec(13), E_vec(12));
for m = N_AlSb+2:N_AlSb+N_InAs+1
    [lattice, E_vec, tao] = parametergen_Hca_InAs(m);
    matrix( 20*m-9:20*m, 20*m+1:20*m+10 ) = Hca(1, 1, lattice(5), lattice(2), lattice(3), lattice(4), E_vec(9), E_vec(10), E_vec(11), E_vec(14), E_vec(15), E_vec(13), E_vec(12)); 
end

m = N_AlSb + N_InAs + 2;
[lattice, E_vec, tao] = parametergen_Hca_IF2(m);
matrix( 1:10, 20*m-9:20*m ) = Hca(1, 1, lattice(5), lattice(2), lattice(3), lattice(4), E_vec(9), E_vec(10), E_vec(11), E_vec(14), E_vec(15), E_vec(13), E_vec(12))'; 

matrix = matrix + matrix';

m = 1;
[lattice1, E_vec1, tao1] = parametergen_Hac_AlSb(m);
[lattice2, E_vec2, tao2] = parametergen_Hca_IF2(N_AlSb + N_InAs + 2);
E_vec = (E_vec1 + E_vec2)/2;
matrix( 20*m-19:20*m-10, 20*m-19:20*m-10 ) = H_diag(E_vec(1), E_vec(3), E_vec(7), E_vec(5));
for m = 2:N_AlSb
   [lattice1, E_vec1, tao1] = parametergen_Hac_AlSb(m);
   [lattice2, E_vec2, tao2] = parametergen_Hca_AlSb(m-1);
   E_vec = (E_vec1 + E_vec2)/2;
   matrix( 20*m-19:20*m-10, 20*m-19:20*m-10 ) = H_diag(E_vec(1), E_vec(3), E_vec(7), E_vec(5));
end
m = N_AlSb + 1;
[lattice1, E_vec1, tao1] = parametergen_Hac_IF1(m);
[lattice2, E_vec2, tao2] = parametergen_Hca_AlSb(m-1);
E_vec = (E_vec1 + E_vec2)/2;
matrix( 20*m-19:20*m-10, 20*m-19:20*m-10 ) = H_diag(E_vec(1), E_vec(3), E_vec(7), E_vec(5));
m = N_AlSb + 2;
[lattice1, E_vec1, tao1] = parametergen_Hac_InAs(m);
[lattice2, E_vec2, tao2] = parametergen_Hca_IF1(m-1);
E_vec = (E_vec1 + E_vec2)/2;
matrix( 20*m-19:20*m-10, 20*m-19:20*m-10 ) = H_diag(E_vec(1), E_vec(3), E_vec(7), E_vec(5));
for m = N_AlSb+3:N_AlSb+N_InAs+1
   [lattice, E_vec, tao] = parametergen_Hac_InAs(m);
   matrix( 20*m-19:20*m-10, 20*m-19:20*m-10 ) = H_diag(E_vec(1), E_vec(3), E_vec(7), E_vec(5));
end
m = N_AlSb + N_InAs + 2;
[lattice1, E_vec1, tao1] = parametergen_Hac_IF2(m);
[lattice2, E_vec2, tao2] = parametergen_Hca_InAs(m-1);
E_vec = (E_vec1 + E_vec2)/2;
matrix( 20*m-19:20*m-10, 20*m-19:20*m-10 ) = H_diag(E_vec(1), E_vec(3), E_vec(7), E_vec(5));

for m = 1:N_AlSb
   [lattice1, E_vec1, tao1] = parametergen_Hac_AlSb(m);
   [lattice2, E_vec2, tao2] = parametergen_Hca_AlSb(m);
   E_vec = (E_vec1 + E_vec2)/2;
   matrix( 20*m-9:20*m, 20*m-9:20*m ) = H_diag(E_vec(2), E_vec(4), E_vec(8), E_vec(6));
end
m = N_AlSb + 1;
[lattice1, E_vec1, tao1] = parametergen_Hac_IF1(m);
[lattice2, E_vec2, tao2] = parametergen_Hca_IF1(m);
E_vec = (E_vec1 + E_vec2)/2;
matrix( 20*m-9:20*m, 20*m-9:20*m ) = H_diag(E_vec(2), E_vec(4), E_vec(8), E_vec(6));
for m = N_AlSb+2:N_AlSb+N_InAs+1
   [lattice, E_vec, tao] = parametergen_Hac_InAs(m);
   matrix( 20*m-9:20*m, 20*m-9:20*m ) = H_diag(E_vec(2), E_vec(4), E_vec(8), E_vec(6));
end
m = N_AlSb + N_InAs + 2;
[lattice1, E_vec1, tao1] = parametergen_Hac_IF2(m);
[lattice2, E_vec2, tao2] = parametergen_Hca_IF2(m);
E_vec = (E_vec1 + E_vec2)/2;
matrix( 20*m-9:20*m, 20*m-9:20*m ) = H_diag(E_vec(2), E_vec(4), E_vec(8), E_vec(6));

