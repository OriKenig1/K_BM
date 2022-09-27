function w = calc_w(x,h)

[N,Dim] = size(x);
Kh_0 = mvnpdf(zeros(1,Dim), zeros(1,Dim), h^2*eye(Dim));

w = zeros(N,1);

% Calculate g_tilde 

for n = 1:N

    Kh = mvnpdf(x(n,:)-x, zeros(1,Dim), h^2*eye(Dim));
    g_hat = sum(Kh) / N;
    g_tilde_curr = g_hat - ( Kh_0 / N );

    w(n) = g_tilde_curr;

end

w = w ./ sum(w);
