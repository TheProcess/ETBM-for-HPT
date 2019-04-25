function mtx = Hca(g1, g2, beta, cosx, cosy, cosz, Esasc, Esaxc, Exasc, Exaxc, Exayc, Exassc, Essaxc)

% Written by Yajun Wei c2004

s3 = sqrt(3);

mtx_sub = zeros(5);
mtx_sub(1,1) = g1 * Esasc;
mtx_sub(1,2) = -g2 * s3 * cosx * Exasc;
mtx_sub(1,3) = g2 * s3 * cosy * Exasc;
mtx_sub(1,4) = g1 * s3 * cosz * Exasc;
mtx_sub(2,1) = g2 * s3 * cosx * Esaxc;
mtx_sub(2,2) = g1 * ( Exaxc + ( 3 * cosx^2 -1 ) * Exayc );
mtx_sub(2,3) = -g1 * 3 * cosx * cosy * Exayc;
mtx_sub(2,4) = -g2 * 3 * cosx * cosz * Exayc;
mtx_sub(2,5) = g2 * s3 * cosx * Essaxc;
mtx_sub(3,1) = -g2 * s3 * cosy * Esaxc;
mtx_sub(3,2) = -g1 * 3 * cosy * cosx * Exayc;
mtx_sub(3,3) = g1 * ( Exaxc + ( 3 * cosy^2 -1 ) * Exayc );
mtx_sub(3,4) = g2 * 3 * cosy * cosz * Exayc;
mtx_sub(3,5) = -g2 * s3 * cosy * Essaxc;
mtx_sub(4,1) = -g1 * s3 * cosz * Esaxc;
mtx_sub(4,2) = -g2 * 3 * cosz * cosx * Exayc;
mtx_sub(4,3) = g2 * 3 * cosz * cosy * Exayc;
mtx_sub(4,4) = g1 * ( Exaxc + ( 3 * cosz^2 -1 ) * Exayc );
mtx_sub(4,5) = -g1 * s3 * cosz * Essaxc;
mtx_sub(5,2) = -g2 * s3 * cosx * Exassc;
mtx_sub(5,3) = g2 * s3 * cosy * Exassc;
mtx_sub(5,4) = g1 * s3 * cosz * Exassc;

mtx_sub = mtx_sub * beta;

mtx_zero = zeros(5);

mtx = [ mtx_sub, mtx_zero;
    mtx_zero, mtx_sub ];