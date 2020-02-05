clear all;
format long
iset=0;
gset=0;
rng(92923932)
ave_psi=0;
global_psi=[];
load matlab;
%%%%%%%%%%%%%%%%%%%%%%% input parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=0.005;
max_time=1.5;
dx=0.5;
fskip=10;
nu=1;
L0=5;
r=0.683;
e0=1;
ndim=3;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for k=[1:ndim]
%     for i=[1:NY]
%         for j=[1:NX]
%             [psi(k,i,j),iset,gset]=gasdev2([ndim,NY,NX],k,i,j,iset,gset);
%             psi(k,i,j)=psi(k,i,j)*1;
%             ave_psi=ave_psi+psi(k,i,j);
%         end
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for k=[1:ndim]
%     for i=[1:NY]
%         for j=[1:NX]
%             [psi(k,i,j),iset,gset]=gasdev2([ndim,NY,NX],k,i,j,iset,gset);
%             psi(k,i,j)=psi(k,i,j)*1;
%             ave_psi=ave_psi+psi(k,i,j);
%         end
%     end
% end

psi=zeros([size(etas,3),size(etas,1),size(etas,2)]);
for k=[1:size(etas,3)]
    for i=[1:size(etas,1)]
        for j=[1:size(etas,2)]
            psi(k,i,j)=etas(i,j,k);
        end
    end
end

global_psi=zeros([(max_time/dt)/fskip,size(psi,2),size(psi,3)]);
global_psi(1,:,:)=scalar_phi(0,psi);    %parameter1 -> 0: just summation ^2, n: summation ^2 except nth array
tmp_psi=zeros(size(psi));
count=1;
for t=[2:max_time/dt]
    for k=[1:size(psi,1)]
        tmp_scalar=scalar_phi(k,psi);
        for i=[1:size(psi,2)]
            for j=[1:size(psi,3)]
                %%%%%%%%%%%%%%%%%%periodic boundary)
                coeff_k=misorientation_read_shock(size(psi,1),k,r,e0);
                tmp_psi(k,i,j)=psi(k,i,j)+dt*L0*(coeff_k*der_2_psi(psi,k,i,j,dx)+(psi(k,i,j)-psi(k,i,j)^3-2*nu*psi(k,i,j)*tmp_scalar(i,j)));
                if tmp_psi(k,i,j)>0.9999
                    tmp_psi(k,i,j)=0.9999;
                elseif tmp_psi(k,i,j)<0.00001
                    tmp_psi(k,i,j)=0.00001;
                end
            end
        end
    end
    if ~rem(t,fskip)
        count=count+1;
        global_psi(count,:,:)=scalar_phi(0,tmp_psi);
    end
    psi=tmp_psi;
end

X=[0:dx:dx*(size(psi,3)-1)];
Y=[0:dx:dx*(size(psi,2)-1)];
for i=[1:(max_time/dt)/fskip]
    Z=reshape(global_psi(i,:),[size(psi,2),size(psi,3)]);
    surf(X,Y,Z,Z)
    caxis([0 1])
    view(2)
    i
    colorbar
    lighting phong;
%     shading interp;
    if i==1
        pause(5);
    end
    pause(1)
end

function [out,iset,gset]=gasdev2(all_size,k,i,j,iset,gset)
out=0;
if k==1
    out=0.01;
    return 
end
if i<all_size(2)/2 & i>all_size(2)/3 & j>all_size(3)/3 & j<all_size(3)/2
    if rand()<0.5
        out=1;
    end
    
end

if i<all_size(2)*0.9 & i>all_size(2)*0.75 & j>all_size(3)*0.75 & j<all_size(3)*0.9
    if k==2
        out=1;
    end
    
end

end

function [psi_der]=der_2_psi(psi,k,i,j,dx)

if (i+1>size(psi,2))
    x1=psi(k,1,j);
else
    x1=psi(k,i+1,j);
end

if(i-1==0)
    x2=psi(k,size(psi,2),j);
else
    x2=psi(k,i-1,j);
end

if (j+1>size(psi,3))
    y1=psi(k,i,1);
else
    y1=psi(k,i,j+1);
end

if(j-1==0)
    y2=psi(k,i,size(psi,3));
else
    y2=psi(k,i,j-1);
end

psi_der=x1+x2+y1+y2-4*psi(k,i,j);

psi_der=psi_der/(dx^2);
if isnan(psi_der)
    [x1,x2,y1,y2,psi(k,i,j)]
    pause
end
end

function [re]=scalar_phi(case_f,psi)
tmp_psi=zeros([size(psi,2),size(psi,3)]);
for i=[1:size(psi,2)]
    for j=[1:size(psi,3)]
        for k=[1:size(psi,1)]
            if k==case_f
                continue
            end
            tmp_psi(i,j)=tmp_psi(i,j)+psi(k,i,j)^2;
        end
    end
end
re=tmp_psi;
end

function [coeff]=misorientation_read_shock(all,k,r,e0)  %parameter1 : Number of order parameter, parameter2:k the array of orderparameter, constant
coeff=e0*sin(2*pi/2.1/all*k)*(1-r*log(sin(2*pi/2.1/all*k)));
%coeff=e0*sin(2*pi/4)*(1-r*log(sin(2*pi/4)));
coeff=coeff^2;
end