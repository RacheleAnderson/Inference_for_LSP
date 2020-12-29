
function S2 = fit_q_to_AIP(L,a,b,c,tau,P)
    
    fun = L+a*exp(-(c/2).*((tau-b).^2));   
    S2 = (1/length(tau)).*(sum((fun-P).^2)); 
    
end
