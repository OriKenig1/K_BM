function covarsw = calc_covarsw(x,means,weights,gamma,lambda)

    dist = (x-means').*(sqrt(weights.*gamma/lambda));

    covarsw = dist'*dist;

end

