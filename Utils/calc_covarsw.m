function covarsw = calc_covarsw(x,means,omega,gamma,lambda)

    dist = (x-means').*(sqrt(omega.*gamma/lambda));

    covarsw = dist'*dist;

end

