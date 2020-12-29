
function S2 = frobDist(L,a_q,b_q,c_q,c_r,time,C_est)
    
    % compute sum of square residuals from fitting C_est to C_model
    % NB. Call this function if centre frequency F0 = 0
    
    t = time;
    s = t';
    tau_R = t*ones(1,length(s))-ones(length(t),1)*s;
    tau_Q = (t*ones(1,length(s))+ones(length(t),1)*s)/2;
   
    RT = exp(-(c_r/8).*(tau_R).^2);
    QH = L+a_q*exp(-(c_q/2)*((tau_Q-b_q).^2));
    C_model = (RT.*QH);
    
    S2 = sqrt(sum(sum((C_est-C_model).^2))); % sqrt of sum of squared residuals
    
end
