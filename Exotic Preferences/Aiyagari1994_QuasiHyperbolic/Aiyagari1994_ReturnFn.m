function F=Aiyagari1994_ReturnFn(aprime_val, a_val, z_val,alpha,delta,gamma,r)
% A list of the parameters to be used
% alpha
% delta
% gamma
% r


F=-Inf;
w=(1-alpha)*((r+delta)/alpha)^(alpha/(alpha-1));
c=w*z_val+(1+r)*a_val-aprime_val; 
%c=wz+(1+r)a_t-a_{t+1}
if c>0
    if gamma==1
        F=log(c);
    else
        F=(c^(1-gamma) -1)/(1-gamma);
    end
end

end