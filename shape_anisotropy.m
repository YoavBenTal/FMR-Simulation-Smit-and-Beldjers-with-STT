%%
function [Na,Nb,Nc] = shape_anisotropy(a,b,c,type)
    %the function culculate the demagnetizing factors Na,Nb,Nc were a>=b>=c>0
    %in SI units (meaning Na+Nb+Nc=1) for three diferent cases:
    %1. a=b=c tipe sphere (type = 1)
    %2. a>>b>=c tipe rod (type = 2)
    %3. a>=b>>c tipe disk (type = 3)
    %source: Demagnetizing Factors of the General Ellipsoid, Physcal Review, J. A. Osborn, June 1945
    
    if type == 1
        Na=1/3;
        Nb=1/3;
        Nc=1/3;
        
    elseif type == 2
        Na=(b*c/(a^2))*(log(4*a/(b+c))-1);
        Nb=c/(b+c)-0.5*(b*c/(a^2))*log(4*a/(b+c))+b*c*(3*b+c)/(4*(a^2)*(b+c));
        Nc=b/(b+c)-0.5*(b*c/(a^2))*log(4*a/(b+c))+b*c*(b+3*c)/(4*(a^2)*(b+c));
        
    elseif type == 3
        argument=sqrt(1-(b/a)^2);
        K=elliptic_integral_tipe_F(pi/2,argument);
        E=elliptic_integral_tipe_E(pi/2,argument);
        Na=(c/a)*sqrt(1-argument^2)*((K-E)/(argument^2));
        Nb=(c/a)*(E-(1-argument^2)*K)/((argument^2)*sqrt(1-argument^2));
        Nc=1-(c*E)/(a*sqrt(1-argument^2));
    
    end
        
    
end


%%
function [K] = elliptic_integral_tipe_F(phi,k)
    %source: Wikipedia, Elliptic integral
    elliptic_integral = @(theta) 1./sqrt(1-(k.^2).*(sin(theta).^2));
    K = integral(elliptic_integral,0,phi);
end

%%
function [K] = elliptic_integral_tipe_E(phi,k)
    %source: Wikipedia, Elliptic integral
    elliptic_integral = @(theta) sqrt(1-(k.^2).*(sin(theta).^2));
    K = integral(elliptic_integral,0,phi);
end