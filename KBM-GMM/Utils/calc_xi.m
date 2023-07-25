function xi = calc_xi(x, GMMStruct, h, u_prime)

    M = GMMStruct.Order;

    N = length(x);

    covars_h = calc_covars_h(GMMStruct.Covars,h);

    xi = zeros(M,1);

    for m = 1:M

        xi(m) = sum( mvnpdf(x,GMMStruct.Means(:,m)', covars_h(:,:,m)) ) ;

    end

    xi = xi ./ ( u_prime*N + 1e-20 );

    
end

