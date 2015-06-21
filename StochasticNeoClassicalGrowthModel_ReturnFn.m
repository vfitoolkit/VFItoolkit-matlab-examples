function F=StochasticNeoClassicalGrowthModel_ReturnFn(kprime_val,k_val,z_val, gamma, alpha, delta)

F=-Inf;
c=exp(z_val)*k_val^alpha-(kprime_val-(1-delta)*k_val);
if c>0
    if gamma==1
        F=log(c);
    else
        F=(c^(1-gamma))/(1-gamma);
    end
end

end