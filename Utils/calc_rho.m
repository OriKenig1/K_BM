function rho = calc_rho(omega, gamma, xi, sigma_inv, sigma_h_inv)

    rho = omega*gamma*sigma_inv - sigma_h_inv*xi;

end

