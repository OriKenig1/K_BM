function alpha = calc_alpha(lambda, eta, xi, alpha)

    alpha = ( 1 ./ (1 + xi) ) .* ( lambda + alpha*eta);

    alpha = alpha ./ sum(alpha);

end

