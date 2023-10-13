function [y,t] = singen(w,n)
    for i = 1:n
        y(i) = sin(w*(i-2));
    end
    
    t = -1:(n-2);
end

