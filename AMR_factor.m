function [factor]=AMR_factor(phi0)
%phi0 in rad
%the function is measured on my sample
    p1 =  1.3428;
    p2 = -6.3650;
    p3 =  9.6800;
    p4 = -3.4080;
    p5 = -2.1880;
    p6 =  0.2706;
    
    factor = abs(p1*phi0^5+p2*phi0^4+p3*phi0^3+p4*phi0^2+p5*phi0+p6);
end
