function xi = calc_xi(x, GMMStruct, h, u_prime)

    M = GMMStruct.Order;

    covars_h = calc_covars_h(GMMStruct.Covars,h);

    xi = zeros(M,1);

    for m = 1:M

        xi(m) = sum( mvnpdf(x,GMMStruct.Means(:,m)', covars_h(:,:,m)) ) ;

    end

    xi = xi ./ ( u_prime + 1e-200 );

    
end

