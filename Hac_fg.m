function mtx = Hac_fg(kvector, tao)

% Written by Yajun Wei c2004

f1 = ( exp( i * sum( kvector .* tao(1,:) )) + exp( i * sum( kvector .* tao(2,:) )) )/4;
f2 = ( exp( i * sum( kvector .* tao(1,:) )) - exp( i * sum( kvector .* tao(2,:) )) )/4;

mtx_sub = zeros(5);
mtx_sub(1,1) = f1;
mtx_sub(1,2) = f2;
mtx_sub(1,3) = f2;
mtx_sub(1,4) = f1;
mtx_sub(2,1) = f2;
mtx_sub(2,2) = f1;
mtx_sub(2,3) = f1;
mtx_sub(2,4) = f2;
mtx_sub(2,5) = f2;
mtx_sub(3,1) = f2;
mtx_sub(3,2) = f1;
mtx_sub(3,3) = f1;
mtx_sub(3,4) = f2;
mtx_sub(3,5) = f2;
mtx_sub(4,1) = f1;
mtx_sub(4,2) = f2;
mtx_sub(4,3) = f2;
mtx_sub(4,4) = f1;
mtx_sub(4,5) = f1;
mtx_sub(5,2) = f2;
mtx_sub(5,3) = f2;
mtx_sub(5,4) = f1;

mtx_zero = zeros(5);

mtx = [ mtx_sub, mtx_zero;
    mtx_zero, mtx_sub ];

