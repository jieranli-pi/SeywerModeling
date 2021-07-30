function phi0 = inlet(phi1,phi2,phix,T1,Tx)
    phi0=(T1*phi1/2+T1*phi2/2+Tx*phix)/(T1+Tx);
end