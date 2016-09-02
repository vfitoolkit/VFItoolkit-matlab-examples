function F=FiniteHorzStochConsSavings_ReturnFn(aprime,a,Wz,gamma,r,Wj)

%jj: age (index 1:J; indicates ages 1 to 10)

W=Wj+exp(Wz);

c=(1+r)*a+W-aprime;

F=-Inf; %(l by aprime)

if c>0
    F=(c.^(1-gamma)-1)/(1-gamma);
end

end