function mtx = Hca_fg(kvector,tao)

% Written by Yajun Wei c2004

g1 = ( exp( -i * sum( kvector .* tao(3,:) )) + exp( -i * sum( kvector .* tao(4,:) )) )/4;
g2 = ( exp( -i * sum( kvector .* tao(3,:) )) - exp( -i * sum( kvector .* tao(4,:) )) )/4;

mtx_sub = zeros(5);
mtx_sub(1,1) = g1;
mtx_sub(1,2) = g2;
mtx_sub(1,3) = g2;
mtx_sub(1,4) = g1;
mtx_sub(2,1) = g2;
mtx_sub(2,2) = g1;
mtx_sub(2,3) = g1;
mtx_sub(2,4) = g2;
mtx_sub(2,5) = g2;
mtx_sub(3,1) = g2;
mtx_sub(3,2) = g1;
mtx_sub(3,3) = g1;
mtx_sub(3,4) = g2;
mtx_sub(3,5) = g2;
mtx_sub(4,1) = g1;
mtx_sub(4,2) = g2;
mtx_sub(4,3) = g2;
mtx_sub(4,4) = g1;
mtx_sub(4,5) = g1;
mtx_sub(5,2) = g2;
mtx_sub(5,3) = g2;
mtx_sub(5,4) = g1;

mtx_zero = zeros(5);

mtx = [ mtx_sub, mtx_zero;
    mtx_zero, mtx_sub ];