function F=Hansen1985_ReturnFn(d_val, aprime_val, a_val, z_val, alpha, delta, gamma)

F=-Inf;
C=exp(z_val)*(a_val^alpha)*(d_val^(1-alpha))-(aprime_val-(1-delta)*a_val);
if C>0
    F=log(C)-gamma*d_val;
end

end