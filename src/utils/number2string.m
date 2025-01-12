function str = number2string(n)
    if n < 10
        str = strcat('00', string(n));
    elseif n>=10 && n<100
        str = strcat('0', string(n));
    else
        str = string(n);
    end
end
