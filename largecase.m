%All variables are input between parameter input start and input end
%only change these parameters

close
clear
clc

tic
%%%%%%%%%%%%%%%%%%%%%%%%%parameter input start
phiin=-10;      %potential in
phiout=10;    %potential_out
T_wire=58479532; %conductivity of electrical wire
T_air=1/10000000;
%T_air=1;
T_soil=1/50;
%T_soil=1;
T_water=1;
T_pipe=1/(15*10^(15));

%dimension & position
onegrid=0.5;%5mm=0.5cm
%actual
width=80;%80cm
length=100;%100cm
pipe_center=[30,35];%center at 30cm, 35cm
pipe_inradius=10;
pipe_thickness=0.5;
leakage_start=170;%the degree of leakage clockwise starting
leakage_end=190;
electrode_center=[30,35];
electrode_radius=3.5;
pin_start=[1,85];%starting at depth 0cm, length 85cm
pin_end=[30,86];%ending at depth 30cm, width 4cm
water_depth=30;%water depth inside of pipe =30cm
%%%%%%%%%%%%%%%%parameter input end

%transpose
width=width/onegrid;
length=length/onegrid;
pipe_center=pipe_center/onegrid;
pipe_inradius=pipe_inradius/onegrid;
pipe_thickness=pipe_thickness/onegrid;
pipe_outradius=pipe_inradius+pipe_thickness;
electrode_center=electrode_center/onegrid;
electrode_radius=electrode_radius/onegrid;
water_depth=water_depth/onegrid;
pin_start=[pin_start(1,1),pin_start(1,2)/onegrid];
pin_end=pin_end/onegrid;

%prepare workspace
[S_2DT, S_2Dphi]=diworkspace(width,length,T_soil);
[Coefficient_matrix, bvec]=systemprep(S_2Dphi,phiin,phiout);

%prepare T matrix
%ring
S_2DT=water_air(S_2DT,water_depth,T_water,T_air);
S_2DT=outside_ring(S_2DT,pipe_outradius,pipe_center,T_soil);
S_2DT=leaked_ring(S_2DT,pipe_outradius,pipe_inradius,pipe_center,T_pipe,leakage_start,leakage_end);
%S_2DT=sealring(S_2DT,pipe_inradius,pipe_center,0);
S_2DT=ring(S_2DT,electrode_radius,0,electrode_center,T_wire);%electrode
%pin coordinte
S_2DT=pin(S_2DT,pin_start,pin_end,T_wire);

%preparation Coefficientmatrix
Coefficient_matrix=stdprep(Coefficient_matrix,S_2Dphi,S_2DT);
%pin
Coefficient_matrix = pin_potential(Coefficient_matrix,S_2Dphi,pin_start,[pin_end(1)+1,pin_end(2)+1]);
%circular
Coefficient_matrix= cir_electrode_potential(Coefficient_matrix,S_2Dphi,S_2DT,T_wire,electrode_center,electrode_radius);

%solve
phi_vector= sparse(Coefficient_matrix)\sparse(bvec);
S_2Dphi= direshape(phi_vector,S_2Dphi);
disp('Current=')
Current=abs((S_2Dphi(electrode_center(1),electrode_center(2))-phiin)*T_wire)
disp('calculation')
toc
tic
%visualise T matrix
plot_trans(S_2DT,'Transmissitivity')

%visualise potential matrix
plot_phi(S_2Dphi,'potential')

flow_contour(S_2Dphi,electrode_center,electrode_radius,'flow_contour')
disp('plotting')
toc


function [S_2DT, S_2Dphi]=diworkspace(width,len,basisT)
    S_2DT=ones(width,len)*basisT;
    S_2Dphi=sparse(width+1,len+1);
end

function [Coefficient_matrix, bvec]=systemprep(S_phi,phiin,phiout)
    siz=size(S_phi,1)*size(S_phi,2);   
    Coefficient_matrix= speye(siz+2);%+2 for later introduce of in & out
    bvec=sparse(siz+2,1);
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

function S_2DT=water_air(S_2DT,water_depth,T_water,T_air)
for i=1:(water_depth-1)
    S_2DT(i,:)=T_air;
end
for i=water_depth:size(S_2DT,1)
    S_2DT(i,:)=T_water;
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

function S_2DT=sealring(S_2DT,radius,center,value)
for i = 1:size(S_2DT,1)%len
    for j = 1:size(S_2DT,2)%width
        if ceil(((i-center(1))^2+(j-center(2))^2)-(radius+1)^2)<=0
            if ceil(((i-center(1))^2+(j-center(2))^2)-radius^2)>=0
                S_2DT(i,j)=value;
            end
        end
    end
end
%quadracheck
for i = 1:(size(S_2DT,1)-1)%len
    for j = 1:(size(S_2DT,2)-1)%width
        v_1 = [i,j] - center;
        v_2 = [center(1)-radius,center(2)] - center;
        x1=v_1(1);
        y1=v_1(2);
        x2=v_2(1);
        y2=v_2(2);
        Theta = atan2d(x1*y2-y1*x2,x1*x2+y1*y2);
        if Theta>=0 && Theta<=180
            [S_2DT(i,j+1),S_2DT(i,j),S_2DT(i+1,j+1),S_2DT(i+1,j)]=quadracheck(S_2DT(i,j+1),S_2DT(i,j),S_2DT(i+1,j+1),S_2DT(i+1,j),value);
        else
            [S_2DT(i,j),S_2DT(i,j+1),S_2DT(i+1,j),S_2DT(i+1,j+1)]=quadracheck(S_2DT(i,j),S_2DT(i,j+1),S_2DT(i+1,j),S_2DT(i+1,j+1),value);
        end
    end
end
end

function [A,B,C,D]=quadracheck(A,B,C,D,value)
if A==0
    if B~=0
        if C~=0
            if D==0
                C=value;
            end
        end
    end
end
if A~=0
    if B==0
        if C==0
            if D~=0
                A=value;
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
                CosTheta = max(min(dot(v_1,v_2)/(norm(v_1)*norm(v_2)),1),-1);
                Theta = real(acosd(CosTheta));
                if j<center(2)
                    Theta=360-Theta;
                end
                if Theta>=leak_start && Theta<=leak_end
                else
                    S_2DT(i,j)=value;
                end
            end
        end
    end
end
%quadracheck
for i = 1:(size(S_2DT,1)-1)%len
    for j = 1:(size(S_2DT,2)-1)%width
        v_1 = [i,j] - center;
        v_2 = [center(1)-inradius,center(2)] - center;
        x1=v_1(1);
        y1=v_1(2);
        x2=v_2(1);
        y2=v_2(2);
        Theta = atan2d(x1*y2-y1*x2,x1*x2+y1*y2);
        if Theta>=0 && Theta<=180
            [S_2DT(i,j+1),S_2DT(i,j),S_2DT(i+1,j+1),S_2DT(i+1,j)]=quadracheck(S_2DT(i,j+1),S_2DT(i,j),S_2DT(i+1,j+1),S_2DT(i+1,j),value);
        else
            [S_2DT(i,j),S_2DT(i,j+1),S_2DT(i+1,j),S_2DT(i+1,j+1)]=quadracheck(S_2DT(i,j),S_2DT(i,j+1),S_2DT(i+1,j),S_2DT(i+1,j+1),value);
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

function S_2DT=outside_ring(S_2DT,radius,center,value)
for i = 1:size(S_2DT,1)%len
    for j = 1:size(S_2DT,2)%width
        if ceil(((i-center(1))^2+(j-center(2))^2)-radius^2)>0
            S_2DT(i,j)=value;
        end
    end
end
end

function Coefficient_matrix = cir_electrode_potential(Coefficient_matrix,S_2Dphi,S_2DT,T_wire,center,radius)
%center coordinate using S_2Dphi
list_out=[];
list_in=[];
for i = 1:size(S_2DT,1)%len
    for j = 1:size(S_2DT,2)%width
        if ceil(((i-center(1))^2+(j-center(2))^2)-(radius+1)^2)<=0
            if ceil(((i-center(1))^2+(j-center(2))^2)-radius^2)>0
                list_out=[list_out;[i,j]];
            end
        end
        if ceil(((i-center(1))^2+(j-center(2))^2)-radius^2)<=0
            list_in=[list_in;[i,j]];
        end
    end
end
    list_in_phi=list_in;
    for i = 1:size(list_in,1)%len
        list_in_phi=[list_in_phi;[list_in(i,1)+1,list_in(i,2)]];
        list_in_phi=[[list_in_phi;list_in(i,1),list_in(i,2)+1]];
        list_in_phi=[list_in_phi;[list_in(i,1)+1,list_in(i,2)+1]];
    end
    list_in_phi = unique(list_in_phi,'rows');
    loc=[];
    for i = 1:size(list_in_phi,1)%len
        loc = [loc;(list_in_phi(i,1)-1)*size(S_2Dphi,2)+list_in_phi(i,2)];
    end
    Coefficient_matrix(loc,:)=0;
    Coefficient_matrix(loc,end-1)=T_wire;
    for i = 1:size(list_in_phi,1)%len
        bot=[list_in_phi(i,1)+1,list_in_phi(i,2)];
        top=[list_in_phi(i,1)-1,list_in_phi(i,2)];
        left=[list_in_phi(i,1),list_in_phi(i,2)-1];
        right=[list_in_phi(i,1),list_in_phi(i,2)+1];
        bot_lia = sum(ismember(list_in_phi,bot,'rows'));
        top_lia = sum(ismember(list_in_phi,top,'rows'));
        left_lia = sum(ismember(list_in_phi,left,'rows'));
        right_lia = sum(ismember(list_in_phi,right,'rows'));
        if top_lia==0
            %top
            cola=(list_in_phi(i,1)-2)*size(S_2Dphi,2)+list_in_phi(i,2);
            Ta=0.5*(S_2DT(list_in_phi(i,1)-1,list_in_phi(i,2)-1)+S_2DT(list_in_phi(i,1)-1,list_in_phi(i,2)));
            Coefficient_matrix(loc,cola)=Coefficient_matrix(loc,cola)+Ta;
        end
        if bot_lia==0
            %bot
            cola=list_in_phi(i,1)*size(S_2Dphi,2)+list_in_phi(i,2);
            Ta=0.5*(S_2DT(list_in_phi(i,1),list_in_phi(i,2)-1)+S_2DT(list_in_phi(i,1),list_in_phi(i,2)));
            Coefficient_matrix(loc,cola)=Coefficient_matrix(loc,cola)+Ta;
        end
        if left_lia==0
            %left
            colb=(list_in_phi(i,1)-1)*size(S_2Dphi,2)+list_in_phi(i,2)-1;
            Tb=0.5*(S_2DT(list_in_phi(i,1)-1,list_in_phi(i,2)-1)+S_2DT(list_in_phi(i,1),list_in_phi(i,2)-1));
            Coefficient_matrix(loc,colb)=Coefficient_matrix(loc,colb)+Tb;
        end
        if right_lia==0
            %right
            colb=(list_in_phi(i,1)-1)*size(S_2Dphi,2)+list_in_phi(i,2)+1;
            Tb=0.5*(S_2DT(list_in_phi(i,1),list_in_phi(i,2))+S_2DT(list_in_phi(i,1)-1,list_in_phi(i,2)));
            Coefficient_matrix(loc,colb)=Coefficient_matrix(loc,colb)+Tb;
        end
    end
    for i = 1:size(loc)%len
        Coefficient_matrix(loc(i),loc(i))=-sum(Coefficient_matrix(loc(i),:));
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

function plot_trans(S_2DT,name)
figure
X=1:size(S_2DT,2);
Y=1:size(S_2DT,1);
C = sin(S_2DT)./S_2DT;
mesh(X,Y,S_2DT,C)
title(append(name,'transmissitivity'))
savefig(append(name,'transmissitivity.fig'))
end

function plot_phi(S_2Dphi,name)
figure
X=1:size(S_2Dphi,2);
Y=1:size(S_2Dphi,1);
%C = sin(S_2Dphi)./S_2Dphi;
mesh(X,Y,S_2Dphi)
title(append(name,'Potential'))
savefig(append(name,'Potential.fig'))
end

function flow_contour(S_2Dphi,center,radius,name)
length=size(S_2Dphi,2)-1;
width=size(S_2Dphi,1)-1;
[X,Y] = meshgrid(0:length,0:width);
[DX,DY] = gradient(S_2Dphi ,1);
figure
hold on
contour(X,Y,S_2Dphi)
title(append(name,'flowandcontour'))
savefig(append(name,'flowcontour.fig'))
startx = [];
starty = [];
for i = 1:size(S_2Dphi,1)%len
    for j = 1:size(S_2Dphi,2)%width
        if ceil(((i-center(1))^2+(j-center(2))^2)-radius^2)<=0
            if ceil(((i-center(1))^2+(j-center(2))^2)-(radius-1)^2)>=0
                startx=[startx,j];
                starty=[starty,i];
            end
        end
    end
end
streamline(X,Y,DX,DY,startx,starty)
set(gca, 'YDir','reverse')
axis equal
hold off
end
