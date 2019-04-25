function mtx = H_diag(Es, Ep, Delta, Ess)

% Written by Yajun Wei c2004

mtx_11 = [ Es,        0,        0,       0,     0;
            0,       Ep,    -i*Delta/3,  0,     0;
            0,   i*Delta/3,     Ep,      0,     0;
            0,        0,        0,      Ep,     0;
            0,        0,        0,       0,    Ess ];

mtx_spin = [ 0,        0,        0,       0,       0;
             0,        0,        0,    Delta/3,    0;
             0,        0,        0,  -i*Delta/3,   0;
             0,    -Delta/3,  i*Delta/3,  0,       0;
             0,        0,        0,       0,       0 ];

H_zero = zeros(5);

mtx = [ mtx_11, mtx_spin;
    mtx_spin', conj(mtx_11) ];
