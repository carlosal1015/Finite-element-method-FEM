function z = comb(k,n)
    z = 1;
    for i = 0:(k-1)
      z = z * (n-i)/(k-i)
    end
end