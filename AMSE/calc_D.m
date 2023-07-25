function [D,v] = calc_D(x,h,GMMStruct)   

    N = length(x);   

    % Building D

    [w, sum_g_tild] = calc_w(x,h);

    psi = (N-1)*w;
    
    c = GMMStruct.q - GMMStruct.y;

    d = calc_lil_d(x,h,GMMStruct,c,sum_g_tild);

    z = ( GMMStruct.f_h ./ GMMStruct.u ) .* ( GMMStruct.b - GMMStruct.y );

    v =  c.*psi' + d - z;

    D = v*v' ./ N;