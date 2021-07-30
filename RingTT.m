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

T1=Ringmatrix(T_matrix,100,5,3,[5,5]);
function T1= Ringmatrix(T,value,outradius,inradius,center)
T1=T;
rows=size(T,1);
cols=size(T,2);
for i = 1:rows
    for j = 1:cols
        if (i-center(1))^2+(j-center(2))^2>inradius^2 && (i-center(1))^2+(j-center(2))^2<outradius^2
            T(i,j)=value;
        end
    end
end
end
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

function phi1 = BRmultientry(phi0,i,size,T,T_wire,length,width,starti) %starti is for the BRnode, rectangle
    startrow=ceil(starti/size);
    startcol=starti-(startrow-1)*size;
    startrow=startrow-1;
    phi1=phi0;
    phi1(1,i)=T_wire;
    for a = 1:length
        phi1(1,i)=phi1(1,i)+T(startrow-width+1-1,startcol-a+1-1);
        phi1(1,starti-width*size-a+1)=-T(startrow-width+1-1,startcol-a+1-1);
    end
    for a = 1:width
        phi1(1,i)=phi1(1,i)+T(startrow-a+1-1,startcol-length+1-1);
        phi1(1,starti-(a-1)*size-length)=-T(startrow-a+1-1,startcol-length+1-1);
    end
    phi1(1,1)=-T_wire;
end