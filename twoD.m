close
clear
clc
phiin=10;      %potential in
phiout=-10;    %potential_out
T_wire=100; %conductivity of electrical wire

%prepare workspace
width=120;
length=150;
[S_2DT, S_2Dphi]=diworkspace(width,length,1);
[Coefficient_matrix, bvec]=systemprep(S_2Dphi,phiin,phiout);

%prepare T matrix
%pin coordinte
pin_start=[1,75];
pin_end=[5,75];
S_2DT=pin(S_2DT,pin_start,pin_end,0);
%rec_electrode
rec_start=[72,57];
rec_end=[75,60];
S_2DT=pin(S_2DT,rec_start,rec_end,0);
%ring
S_2DT=ring(S_2DT,39,38,[75,60],0);

%visualise T matrix
figure
surf(S_2DT)
title('Transmissitivity')

%preparation Coefficientmatrix
Coefficient_matrix=stdprep(Coefficient_matrix,S_2Dphi,S_2DT);
%pin
Coefficient_matrix = pin_potential(Coefficient_matrix,S_2Dphi,pin_start,pin_end);
%electrode
Coefficient_matrix = fal_electrode_potential(Coefficient_matrix,S_2Dphi,[72,57],[75,60]);

%solve
phi_vector = Coefficient_matrix\bvec;
phi_matrix = direshape(phi_vector,S_2Dphi);

%visualise potential matrix
figure
surf(phi_matrix)
title('fake (fixed) electrode')

%rectangular entry electrode
Coefficient_matrix_rec = rec_electrode_potential(Coefficient_matrix,S_2Dphi,S_2DT,T_wire,[72,57],[75,60]);

%solve
phi_vector_rec = Coefficient_matrix_rec\bvec;
phi_matrix_rec = direshape(phi_vector_rec,S_2Dphi);


%visualise potential matrix
figure
surf(phi_matrix)
title('fake (fixed) electrode')

figure
surf(phi_matrix_rec)
title('rectangular entry electrode')

figure
[X,Y] = meshgrid(0:length,0:width);
[DX,DY] = gradient(phi_matrix ,1);
quiver(X,Y,DX*5,DY*5,'AutoScale','off')
title('fake electrode with fixed value as potential in')

hold on
contour(X,Y,phi_matrix)
streamline(X,Y,DX,DY)
axis equal
hold off

figure
[X,Y] = meshgrid(0:length,0:width);
[DX,DY] = gradient(phi_matrix ,1);
quiver(X,Y,DX*5,DY*5,'AutoScale','off')
title('rectangular electrode entry')

hold on
contour(X,Y,phi_matrix)
streamline(X,Y,DX,DY)
axis equal
hold off

%current intotal
disp('Current total=')
Total_current=(phiin-phi_matrix_rec(rec_start(1),rec_start(2)))*T_wire

function [S_2DT, S_2Dphi]=diworkspace(width,len,basisT)
    S_2DT=ones(width,len)*basisT;
    S_2Dphi=zeros(width+1,len+1);
end

function [Coefficient_matrix, bvec]=systemprep(S_phi,phiin,phiout)
    siz=size(S_phi,1)*size(S_phi,2);
    Coefficient_matrix= diag(ones(siz+2,1));%+2 for later introduce of in & out
    bvec=zeros(siz+2,1);
    bvec(end)=phiout;
    bvec(end-1)=phiin;
end

function S_2DT=pin(S_2DT,start,eind,value)
for i =start(1):eind(1)
    for j =start(2):eind(2)
        S_2DT(i,j)=value;
    end
end
end

function S_2DT=ring(S_2DT,outradius,inradius,center,value)
for i = 1:size(S_2DT,1)%len
    for j = 1:size(S_2DT,2)%width
        if ceil(((i-center(1))^2+(j-center(2))^2)-outradius^2)<=0
            if ceil(((i-center(1))^2+(j-center(2))^2)-inradius^2)>=0
                S_2DT(i,j)=value;
            end
        end
    end
end
end

function S_2DT=leaked_ring(S_2DT,outradius,inradius,center,value,leak_start,leak_end)
for i = 1:size(S_2DT,1)%len
    for j = 1:size(S_2DT,2)%width
        if ceil(((i-center(1))^2+(j-center(2))^2)-outradius^2)<=0
            if ceil(((i-center(1))^2+(j-center(2))^2)-inradius^2)>=0
                v_1 = [i,j] - center;
                v_2 = [center(1)-inradius,center(2)] - center;
                x1=v_1(1);
                y1=v_1(2);
                x2=v_2(1);
                y2=v_2(2);
                Theta = atan2d(x1*y2-y1*x2,x1*x2+y1*y2);
                if Theta>=leak_start && Theta<=leak_end
                else
                    S_2DT(i,j)=value;
                end
            end
        end
    end
end
end

function Coefficient_matrix=stdprep(Coefficient_matrix,S_2Dphi,S_2DT)
for i =1:(size(S_2Dphi,1)*size(S_2Dphi,2))
    row=ceil(i/size(S_2Dphi,2));
    col=i-(row-1)*size(S_2Dphi,2);
    if i==1
        Coefficient_matrix(i,:)=TLcorner(Coefficient_matrix(i,:),i, size(S_2Dphi,2));
    elseif i==size(S_2Dphi,2)
        Coefficient_matrix(i,:)=TRcorner(Coefficient_matrix(i,:),i, size(S_2Dphi,2));
    elseif i==(size(S_2Dphi,1)*size(S_2Dphi,2))+1-size(S_2Dphi,2)
        Coefficient_matrix(i,:)=BLcorner(Coefficient_matrix(i,:),i, size(S_2Dphi,2));
    elseif i==size(S_2Dphi,1)*size(S_2Dphi,2)
        Coefficient_matrix(i,:)=BRcorner(Coefficient_matrix(i,:),i, size(S_2Dphi,2));
    elseif row==1
        Coefficient_matrix(i,:) = Topside(Coefficient_matrix(i,:),i, size(S_2Dphi,2),S_2DT,row,col);
    elseif row==size(S_2Dphi,1)
        Coefficient_matrix(i,:) = Botside(Coefficient_matrix(i,:),i, size(S_2Dphi,2),S_2DT,row,col);
    elseif col==1
        Coefficient_matrix(i,:) = Leftside(Coefficient_matrix(i,:),i, size(S_2Dphi,2),S_2DT,row,col);
    elseif col==size(S_2Dphi,2)
        Coefficient_matrix(i,:) = Rightside(Coefficient_matrix(i,:),i, size(S_2Dphi,2),S_2DT,row,col);
    else
        Coefficient_matrix(i,:) = internode(Coefficient_matrix(i,:),i,size(S_2Dphi,2),S_2DT,row,col);
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
function Coefficient_matrix = pin_potential(Coefficient_matrix,S_2Dphi,pin_start,pin_end)
for i =1:(size(S_2Dphi,1)*size(S_2Dphi,2))
    row=ceil(i/size(S_2Dphi,2));
    col=i-(row-1)*size(S_2Dphi,2);
    if row>=pin_start(1)
        if row<=pin_end(1)
            if col>=pin_start(2)
                if col<=pin_end(2)
                    Coefficient_matrix(i,:)=0;
                    Coefficient_matrix(i,i)=1;
                    Coefficient_matrix(i,end)=-1;
                end
            end
        end
    end
end
end
function Coefficient_matrix = fal_electrode_potential(Coefficient_matrix,S_2Dphi,pin_start,pin_end)
for i =1:(size(S_2Dphi,1)*size(S_2Dphi,2))
    row=ceil(i/size(S_2Dphi,2));
    col=i-(row-1)*size(S_2Dphi,2);
    if row>=pin_start(1)
        if row<=pin_end(1)
            if col>=pin_start(2)
                if col<=pin_end(2)
                    Coefficient_matrix(i,:)=0;
                    Coefficient_matrix(i,i)=1;
                    Coefficient_matrix(i,end-1)=-1;
                end
            end
        end
    end
end
end

function Coefficient_matrix = rec_electrode_potential(Coefficient_matrix,S_2Dphi,S_2DT,T_wire,rec_start,rec_end)
for i =1:(size(S_2Dphi,1)*size(S_2Dphi,2))
    row=ceil(i/size(S_2Dphi,2));
    col=i-(row-1)*size(S_2Dphi,2);
    if row>=rec_start(1)
        if row<=rec_end(1)
            if col>=rec_start(2)
                if col<=rec_end(2)
                    Coefficient_matrix(i,:)=0;
                    Coefficient_matrix(i,end-1)=-T_wire;
                    for j=((rec_start(1)-2)*size(S_2Dphi,2)+rec_start(2)):((rec_start(1)-2)*size(S_2Dphi,2)+rec_end(2))
                        col_j=j-(rec_start(1)-2)*size(S_2Dphi,2);
                        Coefficient_matrix(i,j)=-0.5*(S_2DT(rec_start(1)-1,col_j-1)+S_2DT(rec_start(1)-1,col_j));
                    end
                    for j=(rec_end(1)*size(S_2Dphi,2)+rec_start(2)):(rec_end(1)*size(S_2Dphi,2)+rec_end(2))
                        col_j=j-rec_end(1)*size(S_2Dphi,2);
                        Coefficient_matrix(i,j)=-0.5*(S_2DT(rec_end(1),col_j-1)+S_2DT(rec_end(1),col_j));
                    end
                    for j=((rec_start(1)-1)*size(S_2Dphi,2)+rec_start(2)-1):size(S_2Dphi,2):((rec_end(1)-1)*size(S_2Dphi,2)+rec_start(2)-1)
                        row_j=(j-rec_start(2)+1)/size(S_2Dphi,2)+1;
                        Coefficient_matrix(i,j)=-0.5*(S_2DT(row_j-1,rec_start(2)-1)+S_2DT(row_j,rec_start(2)-1));
                    end
                    for j=((rec_start(1)-1)*size(S_2Dphi,2)+rec_end(2)+1):size(S_2Dphi,2):((rec_end(1)-1)*size(S_2Dphi,2)+rec_end(2))
                        row_j=(j-rec_end(2)-1)/size(S_2Dphi,2)+1;
                        Coefficient_matrix(i,j)=-0.5*(S_2DT(row_j-1,rec_end(2))+S_2DT(row_j,rec_end(2)));
                    end
                    Coefficient_matrix(i,i)=-sum(Coefficient_matrix(i,:));
                end
            end
        end
    end
end
end

function S_2Dphi=direshape(phi_vector,S_2Dphi)
for i =1:size(S_2Dphi,1)
    for j =1:size(S_2Dphi,2)
        loc=(i-1)*size(S_2Dphi,2)+j;
        S_2Dphi(i,j)=phi_vector(loc);
    end
end
end