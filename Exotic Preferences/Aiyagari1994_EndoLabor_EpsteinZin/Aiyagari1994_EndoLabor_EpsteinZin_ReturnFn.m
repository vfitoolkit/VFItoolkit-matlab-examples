function F=Aiyagari1994_EndoLabor_EpsteinZin_ReturnFn(d_val, aprime_val, a_val, z_val,chi,r,w)
% A list of the parameters to be used
% alpha
% delta
% gamma
% r

F=-Inf;
c=w*d_val*z_val+(1+r)*a_val-aprime_val; 
%c=wlz+(1+r)a_t-a_{t+1}
if c>0 && d_val<1
    F=((c^(1-chi)) *((1-d_val)^chi));
end

end