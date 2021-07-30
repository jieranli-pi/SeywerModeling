function phi0 = internode(phi1,phi2,phi3,phi4,T1,T2,T3,T4)
    phi0=((T1+T4)*phi1+(T1+T2)*phi2+(T2+T3)*phi3+(T3+T4)*phi4)/(2*(T1+T2+T3+T4));
end