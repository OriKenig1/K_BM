function u = calc_u(x,GMMStruct,h)

    covars_h = calc_covars_h(GMMStruct.Covars,h);

    gmm = gmdistribution(GMMStruct.Means',covars_h,GMMStruct.Alpha);
    
    f_tilde = pdf(gmm,x);

    u = sum(f_tilde);

end