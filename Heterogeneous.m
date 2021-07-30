clear
clc
phi_in=10;      %potential in
phi_out=-10;    %potential_out

size=10;
%conductivity matrix
T_matrix = ones(size-1, size-1);
b_vector = zeros(size^2,1);
b_vector (1,1)=phi_in;
b_vector (end)=phi_out;
Coefficient_matrix= diag(ones(size^2,1));

%Coefficient_matrix_treatment
for i = 1:size^2
    row=ceil(i/size);
    col=i-(row-1)*size;
    %First
    if i==1
        Coefficient_matrix(i,:) = fixed(Coefficient_matrix(i,:));
    %End
    elseif i==size^2
        Coefficient_matrix(i,:) = fixed(Coefficient_matrix(i,:));
    %Topright corner
    elseif i==size
        Coefficient_matrix(i,:) = TRcorner(Coefficient_matrix(size,:),i, size);
    %Botleft corner
    elseif i==size^2-size+1
        Coefficient_matrix(i,:) = BLcorner(Coefficient_matrix(i,:),i, size);
    %Topside
    elseif i>1 && i<size
        Coefficient_matrix(i,:)= Topside(Coefficient_matrix(i,:),i, size,T_matrix,row,col);
    %Botside
    elseif i>size^2-size+1 && i<size^2
        Coefficient_matrix(i,:)= Botside(Coefficient_matrix(i,:),i, size,T_matrix,row,col);
    %Leftside
    elseif rem(i-1,size)==0
        Coefficient_matrix(i,:)= Leftside(Coefficient_matrix(i,:),i, size,T_matrix,row,col);
    %Rightside
    elseif rem(i,size)==0
        Coefficient_matrix(i,:)= Rightside(Coefficient_matrix(i,:),i, size,T_matrix,row,col);
    else
        Coefficient_matrix(i,:)= internode(Coefficient_matrix(i,:),i, size,T_matrix,row,col);
    end
end
%solve
phi_vector = Coefficient_matrix\b_vector;
phi_matrix = reshape(phi_vector,size,size);
%meshgrid
[X,Y] = meshgrid(0:size-1,0:size-1);
U = [phi_matrix(:,1:(end-1))-phi_matrix(:,2:end),zeros(size,1)];
V = [(phi_matrix(1:(end-1),:)-phi_matrix(2:end,:)).',zeros(size,1)].';
quiver(X,Y,U,V)

function phi1 = TRcorner(phi0,i, size)
    phi1=phi0;
    phi1(1,i)=1;
    phi1(1,i-1)=-0.5;
    phi1(1,i+size)=-0.5;
end
function phi1 = BLcorner(phi0,i, size)
    phi1=phi0;
    phi1(1,i)=1;
    phi1(1,i+1)=-0.5;
    phi1(1,i-size)=-0.5;
end
function phi1 = Topside(phi0,i, size,T,row,col)
    phi1=phi0;
    phi1(1,i-1)=-T(row,col-1)/(2*(T(row,col-1)+T(row,col)));
    phi1(1,i+1)=-T((row),(col))/(2*(T((row),(col-1))+T((row),(col))));
    phi1(1,i+size)=-0.5;
end
function phi1 = Botside(phi0,i, size,T,row,col)
    phi1=phi0;
    phi1(1,i-1)=-T((row-1),(col-1))/(2*(T((row-1),(col-1))+T((row-1),(col))));
    phi1(1,i+1)=-T((row-1),(col))/(2*(T((row-1),(col-1))+T((row-1),(col))));
    phi1(1,i-size)=-0.5;
end
function phi1 = Leftside(phi0,i, size,T,row,col)
    phi1=phi0;
    phi1(1,i-size)=-T((row-1),(col))/(2*(T((row-1),(col))+T((row),(col))));
    phi1(1,i+1)=-0.5;
    phi1(1,i+size)=-T((row),(col))/(2*(T((row-1),(col))+T((row),(col))));
end
function phi1 = Rightside(phi0,i, size,T,row,col)
    phi1=phi0;
    phi1(1,i-size)=-T((row-1),(col-1))/(2*(T((row-1),(col-1))+T((row),(col-1))));
    phi1(1,i-1)=-0.5;
    phi1(1,i+size)=-T((row),(col-1))/(2*(T((row-1),(col-1))+T((row),(col-1))));
end
function phi1 = internode(phi0,i,size,T,row,col)
    phi1=phi0;
    phi1(1,i)=2*(T((row-1),(col-1))+T((row-1),(col))+T((row),(col-1))+T((row),(col)));
    phi1(1,i-size)=-T((row-1),(col-1))-T((row-1),(col));
    phi1(1,i-1)=-T((row-1),(col-1))-T((row),(col-1));
    phi1(1,i+size)=-T((row),(col))-T((row),(col-1));
    phi1(1,i+1)=-T((row-1),(col))-T((row),(col));
end
function phi1 = fixed(phi0)
    phi1=phi0;
end