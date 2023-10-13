function [y,t] = singen(w,n)
    b = [0 sin(w)];
    a = [1 -2*cos(w) 1];
    
    x = zeros(n);
    x(1) = 1;
    y = filter(b, a, x);
    
    t = 1:n;
end

