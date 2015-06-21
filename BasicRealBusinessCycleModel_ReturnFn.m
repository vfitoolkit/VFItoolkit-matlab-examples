function F=BasicRealBusinessCycleModel_ReturnFn(d_val,aprime_val,a_val,z_val,alpha,delta,theta,tau)

F=-Inf;
c=exp(z_val)*(a_val^alpha)*(d_val^(1-alpha))-(aprime_val-(1-delta)*a_val);
if c>0
    F=(((c^theta)*((1-d_val)^(1-theta)))^(1-tau))/(1-tau);
end

end