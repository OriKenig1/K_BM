function lambda = calc_lambda(gamma, weights)

    lambda = (gamma' * weights) + 1e-10;

end

