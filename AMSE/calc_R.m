function  [R,C,D] = calc_R(all_x,h,GMMStruct)

    N = length(all_x);

    GMMStruct = calculate_derivatives(all_x,GMMStruct);

    C = calc_C(all_x,h,GMMStruct);
    
    [D,v] = calc_D(all_x,h,GMMStruct);
    
    tmp = ( (C*N) \ v );

    R = (tmp) * (tmp');
