function phi0 = side(phi1,phi2,phi3,T1,T2)
    phi0=phi1/2+(T1*phi2+T2*phi3)/(2*(T1+T2));
end