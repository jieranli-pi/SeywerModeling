clear
clc
phi_in=10;      %potential in
phi_out=-10;    %potential_out
T_wire=100; %conductivity of electrical wire

size=10;
%conductivity matrix
%T_matrix
%T_matrix = ones(size-1, size-1);
T_matrix = [ones(size-1, size-1-ceil(size/2)) 2*ones(size-1, size-1-(size-1-ceil(size/2)))];

b_vector = zeros(size^2+size,1);
b_vector (1,1)=phi_in;
b_vector (end)=phi_out;
Coefficient_matrix= diag(ones(size^2+size,1));

%Coefficient_matrix_treatment
for i = 1:size^2+size
    row=ceil(i/size);
    col=i-(row-1)*size;
    row=row-1;%because extra row for the wire
    %First row
    if i<size+1
        Coefficient_matrix(i,:) = fixed(Coefficient_matrix(i,:));
    %entry
    elseif i-size==1
        Coefficient_matrix(i,:) = entry(Coefficient_matrix(i,:),i,size,T_matrix,T_wire,row,col);
    %End
    elseif i-size==size^2
        Coefficient_matrix(i,:) = fixed(Coefficient_matrix(i,:));
    %Topright corner
    elseif i-size==size
        Coefficient_matrix(i,:) = TRcorner(Coefficient_matrix(size,:),i, size);
    %Botleft corner
    elseif i-size==size^2-size+1
        Coefficient_matrix(i,:) = BLcorner(Coefficient_matrix(i,:),i, size);
    %Topside
    elseif i-size>1 && i-size<size
        Coefficient_matrix(i,:)= Topside(Coefficient_matrix(i,:),i, size,T_matrix,row,col);
    %Botside
    elseif i-size>size^2-size+1 && i-size<size^2
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
phi_matrix = reshape(phi_vector(size+1:end,1),size,size).';

%meshgrid
[X,Y] = meshgrid(0:size-1,0:size-1);
[DX,DY] = gradient(phi_matrix ,1);
quiver(X,Y,DX,DY)

hold on
contour(X,Y,phi_matrix)
streamline(X,Y,DX,DY)
axis equal
hold off

surf(X,Y,phi_matrix)

function phi1 = TRcorner(phi0,i, size)
    phi1=phi0;
    phi1(1,i)=1;
    phi1(1,i-1)=-0.5;
    phi1(1,i+size)=-0.5;
end
function phi1 = BRcorner(phi0,i, size)
    phi1=phi0;
    phi1(1,i)=1;
    phi1(1,i-1)=-0.5;
    phi1(1,i-size)=-0.5;
end
function phi1 = BLcorner(phi0,i, size)
    phi1=phi0;
    phi1(1,i)=1;
    phi1(1,i+1)=-0.5;
    phi1(1,i-size)=-0.5;
end
function phi1 = TLcorner(phi0,i, size)
    phi1=phi0;
    phi1(1,i)=1;
    phi1(1,i+1)=-0.5;
    phi1(1,i+size)=-0.5;
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
function phi1 = entry(phi0,i,size,T,T_wire,row,col)
    phi1=phi0;
    phi1(1,i)=T_wire+T(row,col);
    phi1(1,i+1)=-0.5*T(row,col);
    phi1(1,i+size)=-0.5*T(row,col);
    phi1(1,1)=-T_wire;
end

function phi1 = TLmultientry(phi0,i,size,T,T_wire,length,width,starti) %starti is for the TLnode, rectangle
    startrow=ceil(starti/size);
    startcol=starti-(startrow-1)*size;
    startrow=startrow-1;
    phi1=phi0;
    phi1(1,i)=T_wire;
    for a = 1:length
        phi1(1,i)=phi1(1,i)+T(startrow+width-1,startcol+a-1);
        phi1(1,starti+width*size+a-1)=-T(startrow+width-1,startcol+a-1);
    end
    for a = 1:width
        phi1(1,i)=phi1(1,i)+T(startrow+a-1,startcol+length-1);
        phi1(1,starti+(a-1)*size+length)=-T(startrow+a-1,startcol+length-1);
    end
    phi1(1,1)=-T_wire;
end
function phi1 = TRmultientry(phi0,i,size,T,T_wire,length,width,starti) %starti is for the TRnode, rectangle
    startrow=ceil(starti/size);
    startcol=starti-(startrow-1)*size;
    startrow=startrow-1;
    phi1=phi0;
    phi1(1,i)=T_wire;
    for a = 1:length
        phi1(1,i)=phi1(1,i)+T(startrow+width-1,startcol-a+1-1);
        phi1(1,starti+width*size-a+1)=-T(startrow+width-1,startcol-a+1-1);
    end
    for a = 1:width
        phi1(1,i)=phi1(1,i)+T(startrow+a-1,startcol-length+1-1);
        phi1(1,starti+(a-1)*size-length)=-T(startrow+a-1,startcol-length+1-1);
    end
    phi1(1,1)=-T_wire;
end
function phi1 = BLmultientry(phi0,i,size,T,T_wire,length,width,starti) %starti is for the BLnode, rectangle
    startrow=ceil(starti/size);
    startcol=starti-(startrow-1)*size;
    startrow=startrow-1;
    phi1=phi0;
    phi1(1,i)=T_wire;
    for a = 1:length
        phi1(1,i)=phi1(1,i)+T(startrow-width+1-1,startcol+a-1);
        phi1(1,starti-width*size+a-1)=-T(startrow-width+1-1,startcol+a-1);
    end
    for a = 1:width
        phi1(1,i)=phi1(1,i)+T(startrow-a+1-1,startcol+length-1);
        phi1(1,starti-(a-1)*size+length)=-T(startrow-a+1-1,startcol+length-1);
    end
    phi1(1,1)=-T_wire;
end