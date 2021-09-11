function F=Aiyagari1994_EndoLabor_ReturnFn(d_val, aprime_val, a_val, z_val,gamma_c,gamma_l,chi,r,w,nonSeperableUtility)
% A list of the parameters to be used
% alpha
% delta
% gamma
% r

F=-Inf;
c=w*d_val*z_val+(1+r)*a_val-aprime_val; 
%c=wlz+(1+r)a_t-a_{t+1}
if c>0 && d_val<1
    if nonSeperableUtility==0
        if gamma_c==1 && gamma_l==1
            F=log(c)+chi*log(1-d_val);
        elseif gamma_c==1
            F=log(c)+chi*((1-d_val)^(1-gamma_l) -1)/(1-gamma_l);
        elseif gamma_l==1
            F=(c^(1-gamma_c) -1)/(1-gamma_c)+chi*log(1-d_val);
        else
            F=(c^(1-gamma_c) -1)/(1-gamma_c) +chi*((1-d_val)^(1-gamma_l) -1)/(1-gamma_l);
        end
    else % nonSeperableUtility==1
        % gamma_c plays the role of gamma
        F=(((c^(1-chi))*((1-d_val)^chi))^(1-gamma_c))/(1-gamma_c);
    end
end

end