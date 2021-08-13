close;
clear;
clc;

%phi_in and phi_out
phiin=1000;
phiout=-1000;

%set up workspace
[S_3DT, S_3Dphi]=triworkspace(50,40,5,1);
[Coefficient_matrix, bvec]=systemprep(S_3Dphi,phiin,phiout);

%%T_matrix prepare
%S_3DT=ring(S_3DT,6,5,[25,25],0);
S_3DT=ring(S_3DT,0.5,0,[25,25],0);
S_3DT=pin(S_3DT,[25,1,2],[25,4,2],0);

%%Coeffiecient matrix prepare
Coefficient_matrix = corner_prep(Coefficient_matrix,S_3Dphi);
Coefficient_matrix = side_line_prep(Coefficient_matrix,S_3Dphi,S_3DT);
Coefficient_matrix = side_face_prep(Coefficient_matrix,S_3Dphi,S_3DT);
Coefficient_matrix = inner_prep(Coefficient_matrix,S_3Dphi,S_3DT);
Coefficient_matrix = pintrode_prep(Coefficient_matrix,S_3Dphi,[25,1,2],[26,5,3]);
Coefficient_matrix = electrode_prep(Coefficient_matrix,S_3Dphi,[25,25,1],[26,26,6]);

%solve
phi_vector = Coefficient_matrix\bvec;
S_3Dphi=trireshape(phi_vector,S_3Dphi);

%plot cross1
figure
surf(S_3DT(:,:,2))

figure
surf(S_3Dphi(:,:,2))

%plot cross2
figure
surf(reshape(S_3DT(25,:,:),[size(S_3DT(25,:,:),2),size(S_3DT(25,:,:),3)]))

figure
surf(reshape(S_3Dphi(25,:,:),[size(S_3Dphi(25,:,:),2),size(S_3Dphi(25,:,:),3)]))

function [S_3DT, S_3Dphi]=triworkspace(len,width,height,basisT)
    S_2DT=ones(len,width)*basisT;
    S_2Dphi=zeros(len+1,width+1);
    S_3DT=S_2DT;
    S_3Dphi=cat(3,S_2Dphi,S_2Dphi);
    for i = 2:height
       S_3DT=cat(3,S_3DT,S_2DT);
       S_3Dphi=cat(3,S_3Dphi,S_2Dphi);
    end
end

function [Coefficient_matrix, bvec]=systemprep(S_phi,phiin,phiout)
    siz=size(S_phi,1)*size(S_phi,2)*size(S_phi,3);
    Coefficient_matrix= diag(ones(siz+2,1));%+2 for later introduce of in & out
    bvec=zeros(siz+2,1);
    bvec(end)=phiout;
    bvec(end-1)=phiin;
end

function S_3DT=pin(S_3DT,start,eind,value)
for i =start(1):eind(1)
    for j =start(2):eind(2)
        for k=start(3):eind(3)
            S_3DT(i,j,k)=value;
        end
    end
end
end

function S_3DT=ring(S_3DT,outradius,inradius,center,value)
for i = 1:size(S_3DT,1)%len
    for j = 1:size(S_3DT,2)%width
        for k = 1:size(S_3DT,3)%height
            if ((i-center(1))^2+(j-center(2))^2)-outradius^2<=0
                if ((i-center(1))^2+(j-center(2))^2)-inradius^2>=0
                    S_3DT(i,j,k)=value;
                end
            end
        end
    end
end
end

function Coefficient_matrix = corner_prep(Coefficient_matrix,S_3Dphi)
%len, wid, height=1, 1, 1
i=1;
Coefficient_matrix(i,i+1)=-1;
Coefficient_matrix(i,i+size(S_3Dphi,1))=-1;
Coefficient_matrix(i,i+size(S_3Dphi,1)*size(S_3Dphi,2))=-1;
Coefficient_matrix(i,i)=3;
%len, wid, height=1, max, 1
i=size(S_3Dphi,1);
Coefficient_matrix(i,i-1)=-1;
Coefficient_matrix(i,i+size(S_3Dphi,1))=-1;
Coefficient_matrix(i,i+size(S_3Dphi,1)*size(S_3Dphi,2))=-1;
Coefficient_matrix(i,i)=3;
%len, wid, height=max, 1, 1
i=size(S_3Dphi,1)*(size(S_3Dphi,2)-1)+1;
Coefficient_matrix(i,i+1)=-1;
Coefficient_matrix(i,i-size(S_3Dphi,1))=-1;
Coefficient_matrix(i,i+size(S_3Dphi,1)*size(S_3Dphi,2))=-1;
Coefficient_matrix(i,i)=3;
%len, wid, height=max, max, 1
i=size(S_3Dphi,1)*size(S_3Dphi,2);
Coefficient_matrix(i,i-1)=-1;
Coefficient_matrix(i,i-size(S_3Dphi,1))=-1;
Coefficient_matrix(i,i+size(S_3Dphi,1)*size(S_3Dphi,2))=-1;
Coefficient_matrix(i,i)=3;
%len, wid, height=1, 1, max
i=size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)+1;
Coefficient_matrix(i,i+1)=-1;
Coefficient_matrix(i,i+size(S_3Dphi,1))=-1;
Coefficient_matrix(i,i-size(S_3Dphi,1)*size(S_3Dphi,2))=-1;
Coefficient_matrix(i,i)=3;
%len, wid, height=1, max, max
i=size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)+size(S_3Dphi,1);
Coefficient_matrix(i,i-1)=-1;
Coefficient_matrix(i,i+size(S_3Dphi,1))=-1;
Coefficient_matrix(i,i-size(S_3Dphi,1)*size(S_3Dphi,2))=-1;
Coefficient_matrix(i,i)=3;
%len, wid, height=max, 1, max
i=size(S_3Dphi,1)*size(S_3Dphi,2)*size(S_3Dphi,3)-size(S_3Dphi,1)+1;
Coefficient_matrix(i,i+1)=-1;
Coefficient_matrix(i,i-size(S_3Dphi,1))=-1;
Coefficient_matrix(i,i-size(S_3Dphi,1)*size(S_3Dphi,2))=-1;
Coefficient_matrix(i,i)=3;
%len, wid, height=max, max, max
i=size(S_3Dphi,1)*size(S_3Dphi,2)*size(S_3Dphi,3);
Coefficient_matrix(i,i-1)=-1;
Coefficient_matrix(i,i-size(S_3Dphi,1))=-1;
Coefficient_matrix(i,i-size(S_3Dphi,1)*size(S_3Dphi,2))=-1;
Coefficient_matrix(i,i)=3;
end

function Coefficient_matrix = side_line_prep(Coefficient_matrix,S_3Dphi,S_3DT)
%lineI length=:, width=1, height=1
for i=(1+1):(size(S_3Dphi,1)-1)
    Coefficient_matrix(i,i)=2*(S_3DT(i-1,1,1)+S_3DT(i,1,1));
    Coefficient_matrix(i,i-1)=-S_3DT(i-1,1,1);
    Coefficient_matrix(i,i+1)=-S_3DT(i,1,1);
    Coefficient_matrix(i,i+size(S_3Dphi,1))=-0.5*(S_3DT(i-1,1,1)+S_3DT(i,1,1));
    Coefficient_matrix(i,i+size(S_3Dphi,1)*size(S_3Dphi,2))=-0.5*(S_3DT(i-1,1,1)+S_3DT(i,1,1));
end
%lineII length=1, width=:, height=1
for i=(1:(size(S_3Dphi,2)-2))*size(S_3Dphi,1)+1
    Coefficient_matrix(i,i)=2*(S_3DT(1,(i-1)/size(S_3Dphi,1),1)+S_3DT(1,(i-1)/size(S_3Dphi,1)+1,1));
    Coefficient_matrix(i,i-size(S_3Dphi,1))=-S_3DT(1,(i-1)/size(S_3Dphi,1),1);
    Coefficient_matrix(i,i+size(S_3Dphi,1))=-S_3DT(1,(i-1)/size(S_3Dphi,1)+1,1);
    Coefficient_matrix(i,i+1)=-0.5*(S_3DT(1,(i-1)/size(S_3Dphi,1),1)+S_3DT(1,(i-1)/size(S_3Dphi,1)+1,1));
    Coefficient_matrix(i,i+size(S_3Dphi,1)*size(S_3Dphi,2))=-0.5*(S_3DT(1,(i-1)/size(S_3Dphi,1),1)+S_3DT(1,(i-1)/size(S_3Dphi,1)+1,1));
end          
%lineIII length=1, width=1, height=:
for i=(1:(size(S_3Dphi,3)-2))*size(S_3Dphi,1)*size(S_3Dphi,2)+1
    Coefficient_matrix(i,i)=2*(S_3DT(1,1,(i-1)/(size(S_3Dphi,1)*size(S_3Dphi,2)))+S_3DT(1,1,(i-1)/(size(S_3Dphi,1)*size(S_3Dphi,2))+1));
    Coefficient_matrix(i,i-size(S_3Dphi,1)*size(S_3Dphi,2))=-S_3DT(1,1,(i-1)/(size(S_3Dphi,1)*size(S_3Dphi,2)));
    Coefficient_matrix(i,i+size(S_3Dphi,1)*size(S_3Dphi,2))=-S_3DT(1,1,(i-1)/(size(S_3Dphi,1)*size(S_3Dphi,2))+1);
    Coefficient_matrix(i,i+1)=-0.5*(S_3DT(1,1,(i-1)/(size(S_3Dphi,1)*size(S_3Dphi,2)))+S_3DT(1,1,(i-1)/(size(S_3Dphi,1)*size(S_3Dphi,2))+1));
    Coefficient_matrix(i,i+size(S_3Dphi,1))=-0.5*(S_3DT(1,1,(i-1)/(size(S_3Dphi,1)*size(S_3Dphi,2)))+S_3DT(1,1,(i-1)/(size(S_3Dphi,1)*size(S_3Dphi,2))+1));
end    
%lineIV length=:, width=max, height=1
for i=(size(S_3Dphi,1)*(size(S_3Dphi,2)-1)+2):(size(S_3Dphi,1)*size(S_3Dphi,2)-1)
    Coefficient_matrix(i,i)=2*(S_3DT(i-size(S_3Dphi,1)*(size(S_3Dphi,2)-1)-1,end,1)+S_3DT(i-size(S_3Dphi,1)*(size(S_3Dphi,2)-1),end,1));
    Coefficient_matrix(i,i-1)=-S_3DT(i-size(S_3Dphi,1)*(size(S_3Dphi,2)-1)-1,end,1);
    Coefficient_matrix(i,i+1)=-S_3DT(i-size(S_3Dphi,1)*(size(S_3Dphi,2)-1),end,1);
    Coefficient_matrix(i,i+size(S_3Dphi,1))=-0.5*(S_3DT(i-size(S_3Dphi,1)*(size(S_3Dphi,2)-1)-1,end,1)+S_3DT(i-size(S_3Dphi,1)*(size(S_3Dphi,2)-1),end,1));
    Coefficient_matrix(i,i+size(S_3Dphi,1)*size(S_3Dphi,2))=-0.5*(S_3DT(i-size(S_3Dphi,1)*(size(S_3Dphi,2)-1)-1,end,1)+S_3DT(i-size(S_3Dphi,1)*(size(S_3Dphi,2)-1),end,1));
end
%lineV length=max, width=:, height=1
for i=(2:(size(S_3Dphi,2)-1))*size(S_3Dphi,1)
    Coefficient_matrix(i,i)=2*(S_3DT(end,i/size(S_3Dphi,1)-1,1)+S_3DT(end,i/size(S_3Dphi,1),1));
    Coefficient_matrix(i,i-size(S_3Dphi,1))=-S_3DT(end,i/size(S_3Dphi,1)-1,1);
    Coefficient_matrix(i,i+size(S_3Dphi,1))=-S_3DT(end,i/size(S_3Dphi,1),1);
    Coefficient_matrix(i,i-1)=-0.5*(S_3DT(end,i/size(S_3Dphi,1)-1,1)+S_3DT(end,i/size(S_3Dphi,1),1));
    Coefficient_matrix(i,i+size(S_3Dphi,1)*size(S_3Dphi,2))=-0.5*(S_3DT(end,i/size(S_3Dphi,1)-1,1)+S_3DT(end,i/size(S_3Dphi,1),1));
end   
%lineVI length=max, width=max, height=:
for i=(2:(size(S_3Dphi,3)-1))*size(S_3Dphi,1)*size(S_3Dphi,2)
    Coefficient_matrix(i,i)=2*(S_3DT(end,end,i/size(S_3Dphi,1)/size(S_3Dphi,2)-1)+S_3DT(end,end,i/size(S_3Dphi,1)/size(S_3Dphi,2)));
    Coefficient_matrix(i,i-size(S_3Dphi,1)*size(S_3Dphi,2))=-S_3DT(end,end,i/size(S_3Dphi,1)/size(S_3Dphi,2)-1);
    Coefficient_matrix(i,i+size(S_3Dphi,1)*size(S_3Dphi,2))=-S_3DT(end,end,i/size(S_3Dphi,1)/size(S_3Dphi,2));
    Coefficient_matrix(i,i-1)=-0.5*(S_3DT(end,end,i/size(S_3Dphi,1)/size(S_3Dphi,2)-1)+S_3DT(end,end,i/size(S_3Dphi,1)/size(S_3Dphi,2)));
    Coefficient_matrix(i,i-size(S_3Dphi,1))=-0.5*(S_3DT(end,end,i/size(S_3Dphi,1)/size(S_3Dphi,2)-1)+S_3DT(end,end,i/size(S_3Dphi,1)/size(S_3Dphi,2)));
end  
%lineVII length=max, width=1, height=:
for i=(1:(size(S_3Dphi,3)-2))*size(S_3Dphi,1)*size(S_3Dphi,2)+size(S_3Dphi,1)
    Coefficient_matrix(i,i)=2*(S_3DT(end,1,(i-size(S_3Dphi,1))/(size(S_3Dphi,1)*size(S_3Dphi,2)))+S_3DT(end,1,(i-size(S_3Dphi,1))/(size(S_3Dphi,1)*size(S_3Dphi,2))+1));
    Coefficient_matrix(i,i-size(S_3Dphi,1)*size(S_3Dphi,2))=-S_3DT(end,1,(i-size(S_3Dphi,1))/(size(S_3Dphi,1)*size(S_3Dphi,2)));
    Coefficient_matrix(i,i+size(S_3Dphi,1)*size(S_3Dphi,2))=-S_3DT(end,1,(i-size(S_3Dphi,1))/(size(S_3Dphi,1)*size(S_3Dphi,2))+1);
    Coefficient_matrix(i,i-1)=-0.5*(S_3DT(end,1,(i-size(S_3Dphi,1))/(size(S_3Dphi,1)*size(S_3Dphi,2)))+S_3DT(end,1,(i-size(S_3Dphi,1))/(size(S_3Dphi,1)*size(S_3Dphi,2))+1));
    Coefficient_matrix(i,i+size(S_3Dphi,1))=-0.5*(S_3DT(end,1,(i-size(S_3Dphi,1))/(size(S_3Dphi,1)*size(S_3Dphi,2)))+S_3DT(end,1,(i-size(S_3Dphi,1))/(size(S_3Dphi,1)*size(S_3Dphi,2))+1));
end  
%lineVIII length=:, width=1, height=max
for i=(2:(size(S_3Dphi,1)-1))+size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)
    Coefficient_matrix(i,i)=2*(S_3DT(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)-1,1,end)+S_3DT(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1),1,end));
    Coefficient_matrix(i,i-1)=-S_3DT(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)-1,1,end);
    Coefficient_matrix(i,i+1)=-S_3DT(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1),1,end);
    Coefficient_matrix(i,i+size(S_3Dphi,1))=-0.5*(S_3DT(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)-1,1,end)+S_3DT(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1),1,end));
    Coefficient_matrix(i,i-size(S_3Dphi,1)*size(S_3Dphi,2))=-0.5*(S_3DT(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)-1,1,end)+S_3DT(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1),1,end));
end
%lineIX length=max, width=:, height=max
for i=size(S_3Dphi,1)*size(S_3Dphi,2)*size(S_3Dphi,3)-size(S_3Dphi,1)*(1:(size(S_3Dphi,2)-2))
    Coefficient_matrix(i,i)=2*(S_3DT(end,(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1))/size(S_3Dphi,1)-1,end)+S_3DT(end,(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1))/size(S_3Dphi,1),end));
    Coefficient_matrix(i,i-size(S_3Dphi,1))=-S_3DT(end,(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1))/size(S_3Dphi,1)-1,end);
    Coefficient_matrix(i,i+size(S_3Dphi,1))=-S_3DT(end,(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1))/size(S_3Dphi,1),end);
    Coefficient_matrix(i,i-1)=-0.5*(S_3DT(end,(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1))/size(S_3Dphi,1)-1,end)+S_3DT(end,(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1))/size(S_3Dphi,1),end));
    Coefficient_matrix(i,i-size(S_3Dphi,1)*size(S_3Dphi,2))=-0.5*(S_3DT(end,(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1))/size(S_3Dphi,1)-1,end)+S_3DT(end,(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1))/size(S_3Dphi,1),end));
end 
%lineX length=1, width=:, height=max
for i=size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)+1+(1:(size(S_3Dphi,2)-2))*size(S_3Dphi,1)
    Coefficient_matrix(i,i)=2*(S_3DT(1,(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)-1)/size(S_3Dphi,1),end)+S_3DT(1,(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)-1)/size(S_3Dphi,1)+1,end));
    Coefficient_matrix(i,i-size(S_3Dphi,1))=-S_3DT(1,(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)-1)/size(S_3Dphi,1),end);
    Coefficient_matrix(i,i+size(S_3Dphi,1))=-S_3DT(1,(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)-1)/size(S_3Dphi,1)+1,end);
    Coefficient_matrix(i,i+1)=-0.5*(S_3DT(1,(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)-1)/size(S_3Dphi,1),end)+S_3DT(1,(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)-1)/size(S_3Dphi,1)+1,end));
    Coefficient_matrix(i,i-size(S_3Dphi,1)*size(S_3Dphi,2))=-0.5*(S_3DT(1,(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)-1)/size(S_3Dphi,1),end)+S_3DT(1,(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)-1)/size(S_3Dphi,1)+1,end));
end 
%lineXI length=1, width=max, height=:
for i=size(S_3Dphi,1)*size(S_3Dphi,2)-size(S_3Dphi,1)+1+size(S_3Dphi,1)*size(S_3Dphi,2)*(1:(size(S_3Dphi,3)-2))
    Coefficient_matrix(i,i)=2*(S_3DT(1,end,(i-1+size(S_3Dphi,1))/size(S_3Dphi,1)/size(S_3Dphi,2)-1)+S_3DT(1,end,(i-1+size(S_3Dphi,1))/size(S_3Dphi,1)/size(S_3Dphi,2)));
    Coefficient_matrix(i,i-size(S_3Dphi,1)*size(S_3Dphi,2))=-S_3DT(1,end,(i-1+size(S_3Dphi,1))/size(S_3Dphi,1)/size(S_3Dphi,2)-1);
    Coefficient_matrix(i,i+size(S_3Dphi,1)*size(S_3Dphi,2))=-S_3DT(1,end,(i-1+size(S_3Dphi,1))/size(S_3Dphi,1)/size(S_3Dphi,2));
    Coefficient_matrix(i,i+1)=-0.5*(S_3DT(1,end,(i-1+size(S_3Dphi,1))/size(S_3Dphi,1)/size(S_3Dphi,2)-1)+S_3DT(1,end,(i-1+size(S_3Dphi,1))/size(S_3Dphi,1)/size(S_3Dphi,2)));
    Coefficient_matrix(i,i-size(S_3Dphi,1))=-0.5*(S_3DT(1,end,(i-1+size(S_3Dphi,1))/size(S_3Dphi,1)/size(S_3Dphi,2)-1)+S_3DT(1,end,(i-1+size(S_3Dphi,1))/size(S_3Dphi,1)/size(S_3Dphi,2)));
end
%lineXII length=:, width=max, height=max
for i=(2:(size(S_3Dphi,1)-1))+size(S_3Dphi,1)*size(S_3Dphi,2)*size(S_3Dphi,3)-size(S_3Dphi,1)
    Coefficient_matrix(i,i)=2*(S_3DT(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)-size(S_3Dphi,1)*(size(S_3Dphi,2)-1)-1,end,end)+S_3DT(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)-size(S_3Dphi,1)*(size(S_3Dphi,2)-1),end,end));
    Coefficient_matrix(i,i-1)=-S_3DT(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)-size(S_3Dphi,1)*(size(S_3Dphi,2)-1)-1,end,end);
    Coefficient_matrix(i,i+1)=-S_3DT(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)-size(S_3Dphi,1)*(size(S_3Dphi,2)-1),end,end);
    Coefficient_matrix(i,i-size(S_3Dphi,1))=-0.5*(S_3DT(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)-size(S_3Dphi,1)*(size(S_3Dphi,2)-1)-1,end,end)+S_3DT(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)-size(S_3Dphi,1)*(size(S_3Dphi,2)-1),end,end));
    Coefficient_matrix(i,i-size(S_3Dphi,1)*size(S_3Dphi,2))=-0.5*(S_3DT(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)-size(S_3Dphi,1)*(size(S_3Dphi,2)-1)-1,end,end)+S_3DT(i-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)-size(S_3Dphi,1)*(size(S_3Dphi,2)-1),end,end));
end
end

function Coefficient_matrix = side_face_prep(Coefficient_matrix,S_3Dphi,S_3DT)
%height=1
for i=(2:(size(S_3Dphi,1)-1))
    for j=((1:(size(S_3Dphi,2)-2))*size(S_3Dphi,1))
        jk=j/size(S_3Dphi,1);
        Coefficient_matrix(i+j,i+j)=-1.25*(S_3DT(i-1,jk,1)+S_3DT(i-1,jk+1,1)+S_3DT(i,jk,1)+S_3DT(i,jk+1,1));
        Coefficient_matrix(i+j,i+j-1)=0.5*(S_3DT(i-1,jk,1)+S_3DT(i-1,jk+1,1));
        Coefficient_matrix(i+j,i+j+1)=0.5*(S_3DT(i,jk,1)+S_3DT(i,jk+1,1));
        Coefficient_matrix(i+j,i+j-size(S_3Dphi,1))=0.5*(S_3DT(i-1,jk,1)+S_3DT(i,jk,1));
        Coefficient_matrix(i+j,i+j+size(S_3Dphi,1))=0.5*(S_3DT(i-1,jk+1,1)+S_3DT(i,jk+1,1));
        Coefficient_matrix(i+j,i+j+size(S_3Dphi,1)*size(S_3Dphi,2))=0.25*(S_3DT(i-1,jk,1)+S_3DT(i-1,jk+1,1)+S_3DT(i,jk,1)+S_3DT(i,jk+1,1));
    end
end
%height =max
for i=(2:(size(S_3Dphi,1)-1))
    for j=((1:(size(S_3Dphi,2)-2))*size(S_3Dphi,1))+size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)
        jk=(j-size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1))/size(S_3Dphi,1);
        Coefficient_matrix(i+j,i+j)=-1.25*(S_3DT(i-1,jk,end)+S_3DT(i-1,jk+1,end)+S_3DT(i,jk,end)+S_3DT(i,jk+1,end));
        Coefficient_matrix(i+j,i+j-1)=0.5*(S_3DT(i-1,jk,end)+S_3DT(i-1,jk+1,end));
        Coefficient_matrix(i+j,i+j+1)=0.5*(S_3DT(i,jk,end)+S_3DT(i,jk+1,end));
        Coefficient_matrix(i+j,i+j-size(S_3Dphi,1))=0.5*(S_3DT(i-1,jk,end)+S_3DT(i,jk,end));
        Coefficient_matrix(i+j,i+j+size(S_3Dphi,1))=0.5*(S_3DT(i-1,jk+1,end)+S_3DT(i,jk+1,end));
        Coefficient_matrix(i+j,i+j-size(S_3Dphi,1)*size(S_3Dphi,2))=0.25*(S_3DT(i-1,jk,end)+S_3DT(i-1,jk+1,end)+S_3DT(i,jk,end)+S_3DT(i,jk+1,end));
    end
end
%len=1
for i=(1:(size(S_3Dphi,2)-2))*size(S_3Dphi,1)+1
    for j=(1:(size(S_3Dphi,3)-2))*(size(S_3Dphi,1)*size(S_3Dphi,2))
        ik=(i-1)/size(S_3Dphi,1);
        jk=j/(size(S_3Dphi,1)*size(S_3Dphi,2));
        Coefficient_matrix(i+j,i+j)=-1.25*(S_3DT(1,ik,jk)+S_3DT(1,ik,jk+1)+S_3DT(1,ik+1,jk)+S_3DT(1,ik+1,jk+1));
        Coefficient_matrix(i+j,i+j-1)=0.5*(S_3DT(1,ik,jk)+S_3DT(1,ik,jk+1));
        Coefficient_matrix(i+j,i+j+1)=0.5*(S_3DT(1,ik+1,jk)+S_3DT(1,ik+1,jk+1));
        Coefficient_matrix(i+j,i+j-size(S_3Dphi,1)*size(S_3Dphi,2))=0.5*(S_3DT(1,ik,jk)+S_3DT(1,ik+1,jk));
        Coefficient_matrix(i+j,i+j+size(S_3Dphi,1)*size(S_3Dphi,2))=0.5*(S_3DT(1,ik,jk+1)+S_3DT(1,ik+1,jk+1));
        Coefficient_matrix(i+j,i+j+size(S_3Dphi,1))=0.25*(S_3DT(1,ik,jk)+S_3DT(1,ik,jk+1)+S_3DT(1,ik+1,jk)+S_3DT(1,ik+1,jk+1));
    end
end
%len=max
for i=(2:(size(S_3Dphi,2)-1))*size(S_3Dphi,1)
    for j=(size(S_3Dphi,1)*size(S_3Dphi,2))*(1:(size(S_3Dphi,3)-2))
        ik=i/size(S_3Dphi,1);
        jk=j/(size(S_3Dphi,1)*size(S_3Dphi,2));
        Coefficient_matrix(i+j,i+j)=-1.25*(S_3DT(end,ik-1,jk)+S_3DT(end,ik-1,jk+1)+S_3DT(end,ik,jk)+S_3DT(end,ik,jk+1));
        Coefficient_matrix(i+j,i+j-1)=0.5*(S_3DT(end,ik-1,jk)+S_3DT(end,ik-1,jk+1));
        Coefficient_matrix(i+j,i+j+1)=0.5*(S_3DT(end,ik,jk)+S_3DT(end,ik,jk+1));
        Coefficient_matrix(i+j,i+j-size(S_3Dphi,1)*size(S_3Dphi,2))=0.5*(S_3DT(end,ik-1,jk)+S_3DT(end,ik,jk));
        Coefficient_matrix(i+j,i+j+size(S_3Dphi,1)*size(S_3Dphi,2))=0.5*(S_3DT(end,ik-1,jk+1)+S_3DT(end,ik,jk+1));
        Coefficient_matrix(i+j,i+j-size(S_3Dphi,1))=0.25*(S_3DT(end,ik-1,jk)+S_3DT(end,ik-1,jk+1)+S_3DT(end,ik,jk)+S_3DT(end,ik,jk+1));
    end
end
%wid=1
for i=(2:(size(S_3Dphi,1)-1))
    for j=(size(S_3Dphi,1)*size(S_3Dphi,2))*(1:(size(S_3Dphi,3)-2))
        jk=j/(size(S_3Dphi,1)*size(S_3Dphi,2));
        Coefficient_matrix(i+j,i+j)=-1.25*(S_3DT(i-1,1,jk)+S_3DT(i-1,1,jk+1)+S_3DT(i,1,jk)+S_3DT(i,1,jk+1));
        Coefficient_matrix(i+j,i+j-size(S_3Dphi,1))=0.5*(S_3DT(i-1,1,jk)+S_3DT(i-1,1,jk+1));
        Coefficient_matrix(i+j,i+j+size(S_3Dphi,1))=0.5*(S_3DT(i,1,jk)+S_3DT(i,1,jk+1));
        Coefficient_matrix(i+j,i+j-size(S_3Dphi,1)*size(S_3Dphi,2))=0.5*(S_3DT(i,1,jk)+S_3DT(i-1,1,jk));
        Coefficient_matrix(i+j,i+j+size(S_3Dphi,1)*size(S_3Dphi,2))=0.5*(S_3DT(i,1,jk+1)+S_3DT(i-1,1,jk+1));
        Coefficient_matrix(i+j,i+j+1)=0.25*(S_3DT(i-1,1,jk)+S_3DT(i-1,1,jk+1)+S_3DT(i,1,jk)+S_3DT(i,1,jk+1));
    end
end
%wid=max
for i=((size(S_3Dphi,1)*(size(S_3Dphi,2)-1)+2):(size(S_3Dphi,1)*size(S_3Dphi,2)-1))
    for j=(size(S_3Dphi,1)*size(S_3Dphi,2))*(1:(size(S_3Dphi,3)-2))
        ik=i-(size(S_3Dphi,1)*(size(S_3Dphi,2)-1));
        jk=j/(size(S_3Dphi,1)*size(S_3Dphi,2));
        Coefficient_matrix(i+j,i+j)=-1.25*(S_3DT(ik-1,end,jk)+S_3DT(ik-1,end,jk+1)+S_3DT(ik,end,jk)+S_3DT(ik,end,jk+1));
        Coefficient_matrix(i+j,i+j-size(S_3Dphi,1))=0.5*(S_3DT(ik-1,end,jk)+S_3DT(ik-1,end,jk+1));
        Coefficient_matrix(i+j,i+j+size(S_3Dphi,1))=0.5*(S_3DT(ik,end,jk)+S_3DT(ik,end,jk+1));
        Coefficient_matrix(i+j,i+j-size(S_3Dphi,1)*size(S_3Dphi,2))=0.5*(S_3DT(ik-1,end,jk)+S_3DT(ik,end,jk));
        Coefficient_matrix(i+j,i+j+size(S_3Dphi,1)*size(S_3Dphi,2))=0.5*(S_3DT(ik-1,end,jk+1)+S_3DT(ik,end,jk+1));
        Coefficient_matrix(i+j,i+j-1)=0.25*(S_3DT(ik-1,end,jk)+S_3DT(ik-1,end,jk+1)+S_3DT(ik,end,jk)+S_3DT(ik,end,jk+1));
    end
end
end

function Coefficient_matrix = inner_prep(Coefficient_matrix,S_3Dphi,S_3DT)
for i=2:(size(S_3Dphi,1)-1)
    for j=2:(size(S_3Dphi,2)-1)
        for k=2:(size(S_3Dphi,3)-1)
            loc=(k-1)*size(S_3Dphi,1)*size(S_3Dphi,2)+(j-1)*size(S_3Dphi,1)+i;
            Coefficient_matrix(loc,loc)=0.75*(S_3DT(i-1,j-1,k-1)+S_3DT(i-1,j,k-1)+S_3DT(i-1,j-1,k)+S_3DT(i-1,j,k)+S_3DT(i,j-1,k-1)+S_3DT(i,j,k-1)+S_3DT(i,j-1,k)+S_3DT(i,j,k));
            Coefficient_matrix(loc,loc-1)=0.25*(S_3DT(i-1,j-1,k-1)+S_3DT(i-1,j,k-1)+S_3DT(i-1,j-1,k)+S_3DT(i-1,j,k));
            Coefficient_matrix(loc,loc+1)=0.25*(S_3DT(i,j-1,k-1)+S_3DT(i,j,k-1)+S_3DT(i,j-1,k)+S_3DT(i,j,k));
            Coefficient_matrix(loc,loc-size(S_3Dphi,1))=0.25*(S_3DT(i-1,j-1,k-1)+S_3DT(i,j-1,k-1)+S_3DT(i-1,j-1,k)+S_3DT(i,j-1,k));
            Coefficient_matrix(loc,loc+size(S_3Dphi,1))=0.25*(S_3DT(i-1,j,k-1)+S_3DT(i,j,k-1)+S_3DT(i-1,j,k)+S_3DT(i,j,k));
            Coefficient_matrix(loc,loc-size(S_3Dphi,1)*size(S_3Dphi,2))=0.25*(S_3DT(i-1,j,k-1)+S_3DT(i,j,k-1)+S_3DT(i-1,j-1,k-1)+S_3DT(i,j-1,k-1));
            Coefficient_matrix(loc,loc+size(S_3Dphi,1)*size(S_3Dphi,2))=0.25*(S_3DT(i-1,j,k)+S_3DT(i,j,k)+S_3DT(i-1,j-1,k)+S_3DT(i,j-1,k));
        end
    end
end
end

function Coefficient_matrix=pintrode_prep(Coefficient_matrix,S_3Dphi,start,eind)
for i =start(1):eind(1)
    for j =start(2):eind(2)
        for k=start(3):eind(3)
            loc=(k-1)*size(S_3Dphi,1)*size(S_3Dphi,2)+(j-1)*size(S_3Dphi,1)+i;
            Coefficient_matrix(loc,:)=0;
            Coefficient_matrix(loc,loc)=1;
            Coefficient_matrix(loc,end)=-1;
        end
    end
end
end

function S_3Dphi=trireshape(phi_vector,S_3Dphi)
for i =1:size(S_3Dphi,1)
    for j =1:size(S_3Dphi,2)
        for k=1:size(S_3Dphi,3)
            loc=(k-1)*size(S_3Dphi,1)*size(S_3Dphi,2)+(j-1)*size(S_3Dphi,1)+i;
            S_3Dphi(i,j,k)=phi_vector(loc);
        end
    end
end
end

function Coefficient_matrix=electrode_prep(Coefficient_matrix,S_3Dphi,start,eind)
for i =start(1):eind(1)
    for j =start(2):eind(2)
        for k=start(3):eind(3)
            loc=(k-1)*size(S_3Dphi,1)*size(S_3Dphi,2)+(j-1)*size(S_3Dphi,1)+i;
            Coefficient_matrix(loc,:)=0;
            Coefficient_matrix(loc,loc)=1;
            Coefficient_matrix(loc,end-1)=-1;
        end
    end
end
end

function Coefficient_matrix=electrode_prepI(Coefficient_matrix,S_3Dphi,S_3DT,T_wire,start,eind)
loc1=(size(S_3Dphi,1)*(start(2)-1)+start(1))-1;
loc2=(size(S_3Dphi,1)*(start(2)-1)+start(1))-size(S_3Dphi,1);   
loc3=size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)+(size(S_3Dphi,1)*(start(2)-1)+start(1))-1;
loc4=size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)+(size(S_3Dphi,1)*(start(2)-1)+start(1))-size(S_3Dphi,1);
loc5=(size(S_3Dphi,1)*(start(2)-1)+eind(1))+1;
loc6=(size(S_3Dphi,1)*(start(2)-1)+eind(1))-size(S_3Dphi,1);   
loc7=size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)+(size(S_3Dphi,1)*(start(2)-1)+eind(1))+1;
loc8=size(S_3Dphi,1)*size(S_3Dphi,2)*(size(S_3Dphi,3)-1)+(size(S_3Dphi,1)*(start(2)-1)+eind(1))-size(S_3Dphi,1); 
loc9=loc1+1*size(S_3Dphi,1);
loc10=loc2+3*size(S_3Dphi,1);
loc11=loc3+1*size(S_3Dphi,1);
loc12=loc4+3*size(S_3Dphi,1);
loc13=loc5+1*size(S_3Dphi,1);
loc14=loc6+3*size(S_3Dphi,1);
loc15=loc7+1*size(S_3Dphi,1);
loc16=loc8+3*size(S_3Dphi,1);
for i =start(1):eind(1)
    for j =start(2):eind(2)
        for k=start(3):eind(3)
            loc=(k-1)*size(S_3Dphi,1)*size(S_3Dphi,2)+(j-1)*size(S_3Dphi,1)+i;
            Coefficient_matrix(loc,:)=0;
            Coefficient_matrix(loc,loc)=T_wire;
            Coefficient_matrix(loc,end-1)=-T_wire;
            Coefficient_matrix(loc,loc1)=-0.5*(S_3DT(start(1)-1,start(2)-1,1)+S_3DT(start(start(1)-1,start(2),1)));
            Coefficient_matrix(loc,loc2)=-0.5*(S_3DT(start(1)-1,start(2)-1,1)+S_3DT(start(start(1),start(2)-1,1)));
            Coefficient_matrix(loc,loc3)=-0.5*(S_3DT(start(1)-1,start(2)-1,end)+S_3DT(start(start(1)-1,start(2),end)));
            Coefficient_matrix(loc,loc4)=-0.5*(S_3DT(start(1)-1,start(2)-1,end)+S_3DT(start(start(1),start(2)-1,end)));
            Coefficient_matrix(loc,loc5)=-0.5*(S_3DT(eind(1),start(2)-1,1)+S_3DT(start(eind(1),start(2),1)));
            Coefficient_matrix(loc,loc6)=-0.5*(S_3DT(eind(1)-1,start(2)-1,1)+S_3DT(start(eind(1),start(2)-1,1)));
            Coefficient_matrix(loc,loc7)=-0.5*(S_3DT(eind(1),start(2)-1,end)+S_3DT(start(eind(1),start(2),end)));
            Coefficient_matrix(loc,loc8)=-0.5*(S_3DT(eind(1)-1,start(2)-1,end)+S_3DT(start(eind(1),start(2)-1,end)));
            Coefficient_matrix(loc,loc9)=-0.5*(S_3DT(start(1)-1,eind(2)-1,1)+S_3DT(start(start(1)-1,eind(2),1)));
            Coefficient_matrix(loc,loc10)=-0.5*(S_3DT(start(1)-1,eind(2),1)+S_3DT(start(start(1),eind(2),1)));
            Coefficient_matrix(loc,loc11)=-0.5*(S_3DT(start(1)-1,eind(2)-1,end)+S_3DT(start(start(1)-1,eind(2),end)));
            Coefficient_matrix(loc,loc12)=-0.5*(S_3DT(start(1)-1,eind(2),end)+S_3DT(start(start(1),eind(2),end)));
            Coefficient_matrix(loc,loc13)=-0.5*(S_3DT(eind(1),eind(2)-1,1)+S_3DT(start(eind(1),eind(2),1)));
            Coefficient_matrix(loc,loc14)=-0.5*(S_3DT(eind(1)-1,eind(2),1)+S_3DT(start(eind(1),eind(2),1)));
            Coefficient_matrix(loc,loc7)=-0.5*(S_3DT(eind(1),eind(2)-1,end)+S_3DT(start(eind(1),eind(2),end)));
            Coefficient_matrix(loc,loc8)=-0.5*(S_3DT(eind(1)-1,eind(2),end)+S_3DT(start(eind(1),eind(2),end)));
            %for kk=(loc1+size(S_3Dphi,1)*size(S_3Dphi,2)):(size(S_3Dphi,1)*size(S_3Dphi,2)):(loc3+(size(S_3Dphi,1)*size(S_3Dphi,2)))
        end
    end
    end
end
%function Coefficient_matrix=electrode_prepI(Coefficient_matrix,S_3Dphi,S_3DT,radius,center,value)
%    for i = 1:size(S_3Dphi,1)%len
%        for j = 1:size(S_3Dphi,2)%width
%            for k = 1:size(S_3Dphi,3)%height
%                if ((i-center(1))^2+(j-center(2))^2)-radius^2<=0
%                    loc=(k-1)*size(S_3Dphi,1)*size(S_3Dphi,2)+(j-1)*size(S_3Dphi,1)+i;
%                    elseif ((i-center(1))^2+(j-center(2))^2)-(radius+1)^2<=0
%                    loc=(k-1)*size(S_3Dphi,1)*size(S_3Dphi,2)+(j-1)*size(S_3Dphi,1)+i;
%                end
%            end
%        end
%    end
%end
%end