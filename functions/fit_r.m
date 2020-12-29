
function S2 = fit_r(c,tau,r)
    
    fun = exp(-(c/8).*tau.^2);
    S2 = (1/length(tau)).*(sum((fun-r).^2)); 
    
end
