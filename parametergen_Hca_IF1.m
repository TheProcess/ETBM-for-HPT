function [lattice_par, E_vector, tao] = parametergen_Hca_IF1(m)

% Written by Yajun Wei c2004

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

ai = a_InSb;

D = D_InSb;

exx = a_para/ai - 1; eyy = exx; ezz = -D * exx; % The strain

a_perp = ai * ( 1 + ezz );

a_diagonal = sqrt( a_perp^2 + 2 * a_para^2 );

cosx = a_para/a_diagonal; cosy = cosx; cosz = a_perp/a_diagonal;

beta  = 3/((1 + exx)^2 + (1 + eyy)^2 + (1 + ezz)^2);

Esa =   Esa_InSb; % 77K
Esc =   Esc_InSb;
Epa =    Epa_InSb;
Epc =    Epc_InSb;
Essa =    Essa_InSb;
Essc =    Essc_InSb;
Delta_a =    Delta_a_InSb;
Delta_c =    Delta_c_InSb;
Esasc =   Esasc_InSb;
Esaxc =    Esaxc_InSb;
Exasc =    Exasc_InSb;
Essaxc =    Essaxc_InSb;
Exassc =    Exassc_InSb;
Exaxc =    Exaxc_InSb;
Exayc =    Exayc_InSb;

tao(1,:) = [ a_para, a_para, a_perp ]/4; tao(2,:) = [ -a_para, -a_para, a_perp ]/4;
tao(3,:) = [ a_para, -a_para, -a_perp ]/4; tao(4,:) = [ -a_para, a_para, -a_perp ]/4;



lattice_par = [a_perp, cosx, cosy, cosz, beta];

E_vector = [Esa, Esc, Epa, Epc, Essa, Essc, Delta_a, Delta_c, Esasc, Esaxc, Exasc, Essaxc, Exassc, Exaxc, Exayc];