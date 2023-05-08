

%% 
 clc; clear all; close all; 

%%


%%
%%
%% 
%for given soil water potentials hb and
%transpiration rates E, we calculate the leaf water potentials

hb=-[10:100:15010]';%hPa bulk soil water potential
E=(10^-6).*[0:1:2400];%transpiration in cm3 s^-1
Rroot= 0.1189*10^7; % hPa cm^-3 s 

r0= 0.05; % cm root radius in cm
r2= 0.8; %bulk+root+rhizo soil radius in cm

%the potential h is given in [cm] and it is
%negative for unsaturated conditions

% van Genuchten-Mualem parameters
h0= -1/0.12;%this is alfa [cm^-1]
l=0.305;
ths=0.51;
thr=0.017;
k0= 2.1e-5;;%cm/s 
tau=2 ;%corresponding to 2+l*(a+2) with l exp for theta(h) and a tortusoity
 
theta = 0.22;


%radial domai
r=[r0:0.01:r2];r=r(:);
dr=(r2-r0)./(size(r)-1);
dr=dr(1);

%parameters for cavitation; 
k0_x=1/Rroot; %hPa-1 cm^3 s-1;
k0_x=k0_x./30; % cm3 hPa-1 d-1 
h0_x=-18000; % cm 
tau_x= 5;

%--> for given soil and leaf water potential we calculate the transpiration
%rate based on the conductivities. 
hroot=zeros(length(hb),length(E));
hx=hroot;
hx_max=hroot;
hleaf=hroot;
Emax=hroot;
hbmat=hroot;
Emat=hroot;

L = 1000 ;

for j1=1:length(hb)%size of psi_soil
for j2=1:length(E)%size of E
    
            Emat(j1,j2)=E(j2); 
            hbmat(j1,j2)=hb(j1);
            csoil=-2*pi*r0*L* k0/(1-tau)/(h0^(-tau))/(r0/2- r0*r2^2*(log(r2)-log(r0))/(r2^2-r0^2) ); %
    
            hroot(j1,j2)=-abs(-E(j2)/csoil+hb(j1)^(1-tau)).^(1/(1-tau));            
            if (hb(j1)-hroot(j1,j2))>200000 
            break
            end
            
            p1 = 0.61267;
            p2 = -0.081536;
            
            %hroot(j1,j2)=hroot(j1,j2)*p1-p2*10^4;
            
            %dissipation in the root
            hx(j1,j2)=hroot(j1,j2)-Rroot.*E(j2);     
            
            %hmax - this is when the system is linear
            hx_max(j1,j2)=hbmat(j1,j2)-Rroot.*E(j2);   
            
            %dissipation in the xylem including cavitation 
            cx=-k0_x/(1-tau_x)*(h0_x^tau_x);
            hleaf(j1,j2)=-abs(E(j2)/cx+hx(j1,j2)^(1-tau_x)).^(1/(1-tau_x));
            if (hx(j1,j2)-hleaf(j1,j2))>30000 
            break
            end
                     
end
end



[hleaf_reg,hb_reg]=meshgrid((-20000:100:0),hb); % Regular grid
E_reg=griddata(hleaf,hbmat,Emat,hleaf_reg,hb_reg);         % 3D Interpolates
Emax_reg=griddata(hx_max,hbmat,Emat,hleaf_reg,hb_reg);         % 3D Interpolates
diff_E=Emax_reg-E_reg;                                  % Diff. with Emax


[FX,FY]=gradient(E_reg,100,10);%FX=derivative dE/dhleaf (dhleaf is 100 in the regualr grid)
%rel_grad=abs(FX)./max(abs(FX),[],2);%for ach psi soil, we derive the gradient along psileaf by the maximum of the gradient for this psi_soil
dE=griddata(hleaf,hbmat,((Emat)./(abs(hleaf))  ),hleaf_reg,hb_reg);         % 3D Interpolates


hb_reg_uni=flipud(unique(hb_reg));
hleaf_reg_uni=(unique(hleaf_reg));
thrs_old=0.05*max(E_reg(:));% threshold to split between the Couvreur and Gardner
gradmax=0.50 ;%criterion for the relative gradient
trajectory2=zeros(length(hb_reg_uni),3) ; % Trajectory between zones


for i=1:length(hb_reg_uni)
    rel_grad(i,:)=abs(FX(i,:))./max(abs(FX(i,:)),[],2);

    pos_all=min(find(rel_grad(i,:)>=gradmax)); 

    if ~isempty(pos_all)
        pos=pos_all(end);
        trajectory2(i,:)=[hb_reg_uni(i),hleaf_reg_uni(pos),E_reg(i,pos)];
    else
        trajectory2(i,:)=[hb_reg_uni(i),hleaf_reg_uni(1),E_reg(i,1)];
    end
   
    
    target=hleaf_reg(i,:);
    pos=find(target>=hb_reg_uni(i));
    E_reg(i,pos)=NaN; % impossible zone
    
end


    zon_reg=ones(size(E_reg)); % Separate zones
for i=1:length(hb_reg_uni)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define zones
    target=hleaf_reg(i,:);
    pos=(target>=hb_reg_uni(i));
    zon_reg(i,pos)=0; % imposible zone
    
    target=trajectory2(i,2);
    test=hleaf_reg(i,:);
    pos=find(test<target);
    zon_reg(i,pos)=2;
end


trajecthleaf1=[hleaf_reg_uni(:)*0+hb_reg_uni(1),hleaf_reg_uni(:),E_reg(1,:)'];
trajecthleaf10=[hleaf_reg_uni(:)*0+hb_reg_uni(2),hleaf_reg_uni(:),E_reg(2,:)'];
trajecthleaf20=[hleaf_reg_uni(:)*0+hb_reg_uni(13),hleaf_reg_uni(:),E_reg(13,:)'];
trajecthleaf30=[hleaf_reg_uni(:)*0+hb_reg_uni(20),hleaf_reg_uni(:),E_reg(20,:)'];
trajecthleaf40=[hleaf_reg_uni(:)*0+hb_reg_uni(22),hleaf_reg_uni(:),E_reg(22,:)'];
trajecthleaf50=[hleaf_reg_uni(:)*0+hb_reg_uni(40),hleaf_reg_uni(:),E_reg(40,:)'];
trajecthleaf70=[hleaf_reg_uni(:)*0+hb_reg_uni(45),hleaf_reg_uni(:),E_reg(45,:)'];


trajecthsoil1=[hb_reg_uni(:),hb_reg_uni(:)*0+hleaf_reg_uni(1),E_reg(:,1)];
trajecthsoil10=[hb_reg_uni(:),hb_reg_uni(:)*0+hleaf_reg_uni(51),E_reg(:,51)];
trajecthsoil50=[hb_reg_uni(:),hb_reg_uni(:)*0+hleaf_reg_uni(100),E_reg(:,100)];
trajecthsoil100=[hb_reg_uni(:),hb_reg_uni(:)*0+hleaf_reg_uni(151),E_reg(:,151)];
trajecthsoil1000=[hb_reg_uni(:),hb_reg_uni(:)*0+hleaf_reg_uni(180),E_reg(:,180)];
trajecthsoil1500=[hb_reg_uni(:),hb_reg_uni(:)*0+hleaf_reg_uni(200),E_reg(:,200)];

%%
%%
figure(100) 

   index = (0:0.35:3.8500);
myMap = colormap(jet(length(index)));
set(0,'DefaultAxesColorOrder',myMap);

hold on
plot(-trajectory2(:,2).*10^-4,trajectory2(:,3),'color',[0.9019,0.5959,0.3459],'linewidth',2, 'Markersize', 9) % the SOL
plot(-trajecthleaf10(166:201,2).*10^-4,trajecthleaf10(166:201,3),'k','linewidth',1)

xlabel ('-{\it\psi_{leaf-x}} [MPa]'); ylabel ('{\itE} [cm^3 s^-^1]');
box on; 
axis([ 0 2 0 2.5e-3])
set(gca,'linewidth',2,'FontSize',15,'fontweight','bold')
set(gca, 'color', 'none')
%%
%%
Trans_OneRoot_SOL_X = -trajectory2(:,2).*10^-4;
Trans_OneRoot_SOL_Y = trajectory2(:,3);

Trans_OneRoot_X = -trajecthleaf10(166:201,2).*10^-4;
Trans_OneRoot_Y = trajecthleaf10(166:201,3);

%'color',[0.9019,0.5959,0.3459]

%%
L = 5000 ;
 
for j1=1:length(hb)%size of psi_soil
for j2=1:length(E)%size of E
    
            Emat(j1,j2)=E(j2); 
            hbmat(j1,j2)=hb(j1);
            csoil=-2*pi*r0*L* k0/(1-tau)/(h0^(-tau))/(r0/2- r0*r2^2*(log(r2)-log(r0))/(r2^2-r0^2) ); %
    
            hroot(j1,j2)=-abs(-E(j2)/csoil+hb(j1)^(1-tau)).^(1/(1-tau));            
            if (hb(j1)-hroot(j1,j2))>200000 
            break
            end
            
            p1 = 0.61267;
            p2 = -0.081536;
                        
            %dissipation in the root
            hx(j1,j2)=hroot(j1,j2)-Rroot.*E(j2);      % 
            
            %hmax - this is when the system is linear
            hx_max(j1,j2)=hbmat(j1,j2)-Rroot.*E(j2);   
            
            %dissipation in the xylem including cavitation 
            cx=-k0_x/(1-tau_x)*(h0_x^tau_x);
            hleaf(j1,j2)=-abs(E(j2)/cx+hx(j1,j2)^(1-tau_x)).^(1/(1-tau_x));
            if (hx(j1,j2)-hleaf(j1,j2))>30000 
            break
            end
                              
end
end


[hleaf_reg,hb_reg]=meshgrid((-20000:100:0),hb); % Regular grid
E_reg=griddata(hleaf,hbmat,Emat,hleaf_reg,hb_reg);         % 3D Interpolates
Emax_reg=griddata(hx_max,hbmat,Emat,hleaf_reg,hb_reg);         % 3D Interpolates
diff_E=Emax_reg-E_reg;                                  % Diff. with Emax

[FX,FY]=gradient(E_reg,100,10);%FX=derivative dE/dhleaf (dhleaf is 100 in the regualr grid)
dE=griddata(hleaf,hbmat,((Emat)./(abs(hleaf))  ),hleaf_reg,hb_reg);         % 3D Interpolates

hb_reg_uni=flipud(unique(hb_reg));
hleaf_reg_uni=(unique(hleaf_reg));
thrs_old=0.05*max(E_reg(:));% threshold to split between the Couvreur and Gardner
gradmax=0.50 ;%criterion for the relative gradient
trajectory2=zeros(length(hb_reg_uni),3) ; % Trajectory between zones

for i=1:length(hb_reg_uni)
    rel_grad(i,:)=abs(FX(i,:))./max(abs(FX(i,:)),[],2);
    pos_all=min(find(rel_grad(i,:)>=gradmax)); 

    if ~isempty(pos_all)
        pos=pos_all(end);
        trajectory2(i,:)=[hb_reg_uni(i),hleaf_reg_uni(pos),E_reg(i,pos)];
    else
        trajectory2(i,:)=[hb_reg_uni(i),hleaf_reg_uni(1),E_reg(i,1)];
    end
      
    target=hleaf_reg(i,:);
    pos=find(target>=hb_reg_uni(i));
    E_reg(i,pos)=NaN; % impossible zone
    
end

    zon_reg=ones(size(E_reg)); % Separate zones
for i=1:length(hb_reg_uni)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    target=hleaf_reg(i,:);
    pos=(target>=hb_reg_uni(i));
    zon_reg(i,pos)=0; % 
    
    target=trajectory2(i,2);
    test=hleaf_reg(i,:);
    pos=find(test<target);
    zon_reg(i,pos)=2;
end


trajecthleaf1=[hleaf_reg_uni(:)*0+hb_reg_uni(1),hleaf_reg_uni(:),E_reg(1,:)'];
trajecthleaf10=[hleaf_reg_uni(:)*0+hb_reg_uni(5),hleaf_reg_uni(:),E_reg(5,:)'];
trajecthleaf20=[hleaf_reg_uni(:)*0+hb_reg_uni(11),hleaf_reg_uni(:),E_reg(11,:)'];
trajecthleaf30=[hleaf_reg_uni(:)*0+hb_reg_uni(15),hleaf_reg_uni(:),E_reg(15,:)'];
trajecthleaf40=[hleaf_reg_uni(:)*0+hb_reg_uni(22),hleaf_reg_uni(:),E_reg(22,:)'];
trajecthleaf50=[hleaf_reg_uni(:)*0+hb_reg_uni(40),hleaf_reg_uni(:),E_reg(40,:)'];
trajecthleaf70=[hleaf_reg_uni(:)*0+hb_reg_uni(45),hleaf_reg_uni(:),E_reg(45,:)'];


trajecthsoil1=[hb_reg_uni(:),hb_reg_uni(:)*0+hleaf_reg_uni(1),E_reg(:,1)];
trajecthsoil10=[hb_reg_uni(:),hb_reg_uni(:)*0+hleaf_reg_uni(51),E_reg(:,51)];
trajecthsoil50=[hb_reg_uni(:),hb_reg_uni(:)*0+hleaf_reg_uni(100),E_reg(:,100)];
trajecthsoil100=[hb_reg_uni(:),hb_reg_uni(:)*0+hleaf_reg_uni(151),E_reg(:,151)];
trajecthsoil1000=[hb_reg_uni(:),hb_reg_uni(:)*0+hleaf_reg_uni(180),E_reg(:,180)];
trajecthsoil1500=[hb_reg_uni(:),hb_reg_uni(:)*0+hleaf_reg_uni(200),E_reg(:,200)];

%%
%%
figure(100)
hold on

 index = (0:0.35:3.8500);
myMap = colormap(jet(length(index)));
set(0,'DefaultAxesColorOrder',myMap);

hold on
plot(-trajectory2(:,2).*10^-4,trajectory2(:,3),'-','color', [0.0909,0.3259,0.3849],'linewidth',2, 'Markersize', 9) % the SOL
plot(-trajecthleaf1(:,2).*10^-4,trajecthleaf1(:,3),'--k','linewidth',1)

 hold on
 set(gca, 'xtick', [ 0 0.5 1 1.5 2])

 hold on
set(gcf,'unit','inch','position',[1.5,1.5, 6,4])
set(gca,'linewidth',2,'FontSize',15,'fontweight','bold')
%%
Trans_FiveRoot_SOL_X = -trajectory2(:,2).*10^-4;
Trans_FiveRoot_SOL_Y = trajectory2(:,3);

Trans_FiveRoot_X = -trajecthleaf10(166:201,2).*10^-4;
Trans_FiveRoot_Y = trajecthleaf10(166:201,3);
