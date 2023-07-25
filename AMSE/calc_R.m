function  [R,C,D] = calc_R(all_x,h,GMMStruct)

    N = length(all_x);

    GMMStruct = calculate_derivatives(all_x,GMMStruct);

    C = calc_C(all_x,h,GMMStruct);
    
    [D,v] = calc_D(all_x,h,GMMStruct);
    
    %{
    cond_C = cond(C);
    if(cond_C > 30)
        %warning("S with large condition number")
        C = C + 0.01 * 10^floor(log10(trace(C))) * log10(cond_C) * eye(GMMStruct.len);
    end

    C = (C+C')./2;
    %}
    
    tmp = ( (C*N) \ v );

    R = (tmp) * (tmp');