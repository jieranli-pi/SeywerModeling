clear
clc
phi_in=10;      %potential in
phi_out=-10;    %potential_out
T_wire=1000000; %conductivity of electrical wire

size=20;
T_matrix = ones(size-1);
phi_matrix = zeros(size);
phi_matrix (end)=phi_out;
diff_count=100; %to count the difference of each time iteration

while diff_count > 0.000001 
    A=phi_matrix(2);
    for i = 1:numel(phi_matrix)
        row=ceil(i/size);
        col=i-(row-1)*size;
        if row == 1 && col ==1
            phi_matrix(i)=inlet(phi_matrix(row,col+1),phi_matrix(row+1,col),phi_in,T_matrix(i),T_wire);
        elseif row ==1 && col ~= 1 && col ~= size
            phi_matrix(i)=side(phi_matrix(row+1,col),phi_matrix(row,col-1),phi_matrix(row,col+1),T_matrix(1,(col-1)),T_matrix(1,col));
        elseif row ==1 && col ==size
            phi_matrix(i)=cornerp(phi_matrix(row,col-1),phi_matrix(row+1,col));
        elseif col ==1 && row ~= 1 && row ~= size
            phi_matrix(i)=side(phi_matrix(row,col+1),phi_matrix(row-1,col),phi_matrix(row+1,col),T_matrix((row-1),1),T_matrix(row,1));
        elseif col ==1 && row == size
            phi_matrix(i)=cornerp(phi_matrix(row-1,col),phi_matrix(row,col+1));
        elseif col ==size && row ~= 1 && row ~= size
            phi_matrix(i)=side(phi_matrix(row,col-1),phi_matrix(row-1,col),phi_matrix(row+1,col),T_matrix((row-1),(size-1)),T_matrix(row,(size-1)));
        elseif row ==size && col ~= 1 && col ~= size
            phi_matrix(i)=side(phi_matrix(row-1,col),phi_matrix(row,col-1),phi_matrix(row,col+1),T_matrix((size-1),(col-1)),T_matrix((size-1),col));
        elseif row ==size && col==size
            phi_matrix(i)=phi_out;
        else
            phi_matrix(i)=internode(phi_matrix(row,col-1),phi_matrix(row,col+1),phi_matrix(row-1,col),phi_matrix(row+1,col),T_matrix((row-1),(col-1)),T_matrix((row-1),(col)),T_matrix((row),(col-1)),T_matrix((row),(col)));
        end
    end
    B=phi_matrix(2);
    diff_count=abs(A-B);
end


function phi0 = side(phi1,phi2,phi3,T1,T2)
    phi0=phi1/2+(T1*phi2+T2*phi3)/(2*(T1+T2));
end

function phi0 = cornerp(phi1,phi2)
    phi0=(phi1+phi2)/2;
end

function phi0 = internode(phi1,phi2,phi3,phi4,T1,T2,T3,T4)
    phi0=((T1+T4)*phi1+(T1+T2)*phi2+(T2+T3)*phi3+(T3+T4)*phi4)/(2*(T1+T2+T3+T4));
end

function phi0 = inlet(phi1,phi2,phix,T1,Tx)
    phi0=(T1*phi1/2+T1*phi2/2+Tx*phix)/(T1+Tx);
end