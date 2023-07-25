function covars_phi = calc_covarsphi(x,means,covars_h)

    phi = mvnpdf(x,means', covars_h) + 1e-10;

    dist = (x-means').*sqrt(phi./sum(phi));

    covars_phi = dist'*dist;

end

