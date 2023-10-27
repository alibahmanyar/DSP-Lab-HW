function y = myconv(h,x)
    M = length(h);
    L = length(x);
    N = M + L - 1;
    X = [x, zeros(1,M)];
    H = [h, zeros(1,L)];
    y = zeros(1, N);
    if(M > L)
        temp = X;
        X = H;
        H = temp;
    end
    for i = 1 : N
        for j = 1 : min(L, M)
            if (i - j + 1 > 0)
                y(i) = y(i) + H(j) * X(i - j + 1);
            end
        end
    end
end