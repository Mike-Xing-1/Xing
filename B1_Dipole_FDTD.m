%---------------Dipole antenna simulation----------------------------------
clear all
close all
clc
NX=100;               % whole space x dimension
NY=100;               % whole space y dimension
NZ=100;               % whole space z dimension
dipole_length=0.45;    % dipole length (wavelength)
radi=0.01;            % dipole radius (wavelength)
cell_x=0.5/(34*2);        % cell dimension 
Nd=floor(dipole_length/cell_x)+1;

xc=ceil(NX/2);
yc=ceil(NY/2);
zc=ceil(NZ/2);

c0=3e8;
dt=cell_x/(2*c0);     %Time step
Z0=377;

px=7;                 %PML length in x diection
py=7;                 %PML length in y diection
pz=7;                 %PML length in z diection

pml_length=px;
p0=0.33;

pxb=NX-px-1;
pyb=NY-py-1;
pzb=NZ-pz-1;

ex=zeros(NX,NY,NZ);  % Initializing the E-field vector
ey=zeros(NX,NY,NZ);
ez=zeros(NX,NY,NZ);

dx=zeros(NX,NY,NZ);  % Initializing the D vector
dy=zeros(NX,NY,NZ);
dz=zeros(NX,NY,NZ);

hx=zeros(NX,NY,NZ);  % Initializing the H-field vector
hy=zeros(NX,NY,NZ);
hz=zeros(NX,NY,NZ);

cx=ones(NX,NY,NZ);
cy=ones(NX,NY,NZ);
cz=ones(NX,NY,NZ);

% Initialization I is used to store the cumulative conductivity effect
I_x = zeros(NX, NY, NZ);  % I for E_x
I_y = zeros(NX, NY, NZ);  % I for E_y
I_z = zeros(NX, NY, NZ);  % I for E_z

curr_dxl=zeros(px,NY,NZ);
curr_dxh=zeros(px,NY,NZ);
curr_dyl=zeros(NX,py,NZ);
curr_dyh=zeros(NX,py,NZ);
curr_dzl=zeros(NX,NY,pz);
curr_dzh=zeros(NX,NY,pz);

curr_hxl=zeros(px,NY,NZ);
curr_hxh=zeros(px,NY,NZ);
curr_hyl=zeros(NX,py,NZ);
curr_hyh=zeros(NX,py,NZ);
curr_hzl=zeros(NX,NY,pz);
curr_hzh=zeros(NX,NY,pz);

curr_dipole=zeros(35,2048); % Initializing the dipole current


%----------------Add Multi-layered Spherical Load Parameters----------------%

% Define the relative permittivity (epsilon_r), conductivity (sigma), and thickness (in meters) for each layer
layers = {
    % [Relative permittivity, Conductivity (S/m), Thickness (meters)]
    [43.8, 0.413, 6/(34*2)],     % White matter (thickness 8 cm)
    [72.7, 2.225, 1/(34*2)],   % CSF (thickness 1.4 cm)
    [13.4, 0.0827, 1/(34*2)]   % Skull (thickness 1.2 cm)
};
epsilon_0 = 8.854e-12;  % Permittivity of vacuum, in F/m
distance = 20; % In grid units

% Define the center position of the sphere (in grid units)
load_center = [xc, yc + distance, zc]; % Sphere center offset along the X-axis

% Initialize the outer radius of the current layer (in grid units)
load_radius = 0;

% Define parameters from the innermost to the outermost layer, incrementing the outer radius
for layer = 1:length(layers)
    epsilon_r_sphere = layers{layer}(1);
    sigma_sphere = layers{layer}(2);
    layer_thickness = round(layers{layer}(3) / cell_x);  % Convert thickness from meters to number of grid units

    % Update the outer radius of the current layer
    inner_radius = load_radius;
    load_radius = load_radius + layer_thickness;

    % Calculate the distance to the center of the sphere
    for i = 1:NX
        for j = 1:NY
            for k = 1:NZ
                % Calculate the distance to the center of the sphere
                distance_to_center = sqrt((i - load_center(1))^2 + (j - load_center(2))^2 + (k - load_center(3))^2);
                if distance_to_center > inner_radius && distance_to_center <= load_radius

                % If the point is inside the sphere, update the conductivity and permittivity
                epsilon_r(i,j,k) = epsilon_r_sphere;
                sigma(i,j,k) = sigma_sphere;

                % Update the electric field coefficient only inside the sphere
                cx(i,j,k) = 1 / (epsilon_r(i,j,k) + (dt / epsilon_0) * sigma(i,j,k));
                cy(i,j,k) = 1 / (epsilon_r(i,j,k) + (dt / epsilon_0) * sigma(i,j,k));
                cz(i,j,k) = 1 / (epsilon_r(i,j,k) + (dt / epsilon_0) * sigma(i,j,k));

                end
            end
        end
    end
end


% % Define the relative permittivity and conductivity of the sphere
% epsilon_r_sphere = 43.8;  % Relative permittivity of the sphere
% sigma_sphere = 0.413;     % Conductivity of the sphere (S/m)
% distance = 2 * 2;
% 
% % Define the radius and center position of the sphere
% load_radius = 8 * 2;  % Radius of the sphere (in grid units)
% load_center = [xc, yc + distance + load_radius, zc];  % Sphere center position
% 
% % Initialize the relative permittivity and conductivity for the entire space
% epsilon_r = ones(NX, NY, NZ);  % Default relative permittivity, set to air (1)
% sigma = zeros(NX, NY, NZ);     % Default conductivity, set to 0
% 
% % Assign specific conductivity and permittivity in the sphere region and update electric field coefficients
% for i = 1:NX
%     for j = 1:NY
%         for k = 1:NZ
%             % Calculate the distance to the center of the sphere
%             distance_to_center = sqrt((i - load_center(1))^2 + (j - load_center(2))^2 + (k - load_center(3))^2);
%             if distance_to_center <= load_radius
%                 % If the point is inside the sphere, update the conductivity and permittivity
%                 epsilon_r(i,j,k) = epsilon_r_sphere;
%                 sigma(i,j,k) = sigma_sphere;
% 
%                 % Update the electric field coefficient only inside the sphere
%                 cx(i,j,k) = 1 / (epsilon_r(i,j,k) + (dt / epsilon_0) * sigma(i,j,k));
%                 cy(i,j,k) = 1 / (epsilon_r(i,j,k) + (dt / epsilon_0) * sigma(i,j,k));
%                 cz(i,j,k) = 1 / (epsilon_r(i,j,k) + (dt / epsilon_0) * sigma(i,j,k));
%             end
%         end
%     end
% end

%---------------Specify Dipole----------------------------%

for k=zc-floor(Nd/2):zc+floor(Nd/2)
    cz(xc,yc,k)=0.0;
end

cz(xc,yc,zc)=1.0;


%---------------Specify Dipole----------------------------%

%--------------PML boundary condition---------------------%
mx1=zeros(1,NX);
nx1=zeros(1,NX);
mx2=ones(1,NX);
nx2=ones(1,NX);
mx3=ones(1,NX);
nx3=ones(1,NX);

my1=zeros(1,NY);
ny1=zeros(1,NY);
my2=ones(1,NY);
ny2=ones(1,NY);
my3=ones(1,NY);
ny3=ones(1,NY);

mz1=zeros(1,NZ);
nz1=zeros(1,NZ);
mz2=ones(1,NZ);
nz2=ones(1,NZ);
mz3=ones(1,NZ);
nz3=ones(1,NZ);

pml_length=7;
for i=1:pml_length
    ax=(pml_length-i+1)/pml_length;
    ax1=p0*ax^3;
    nx1(i)=ax1;
    nx1(NX-i+1)=ax1;
    mx2(i)=1/(1+ax1);
    mx2(NX-i+1)=1/(1+ax1);
    mx3(i)=(1-ax1)/(1+ax1);
    mx3(NX-i+1)=(1-ax1)/(1+ax1);
    ax=(pml_length-i+0.5)/pml_length;
    ax1=p0*ax^3;
    mx1(i)=ax1;
    mx1(NX-i)=ax1;
    nx2(i)=1/(1+ax1);
    nx2(NX-i)=1/(1+ax1);
    nx3(i)=(1-ax1)/(1+ax1);
    nx3(NX-i)=(1-ax1)/(1+ax1);
end

for j=1:pml_length
    ax=(pml_length-j+1)/pml_length;
    ax1=p0*ax^3;
    ny1(j)=ax1;
    ny1(NX-j+1)=ax1;
    my2(j)=1/(1+ax1);
    my2(NX-j+1)=1/(1+ax1);
    my3(j)=(1-ax1)/(1+ax1);
    my3(NX-j+1)=(1-ax1)/(1+ax1);
    ax=(pml_length-j+0.5)/pml_length;
    ax1=p0*ax^3;
    my1(j)=ax1;
    my1(NX-j)=ax1;
    ny2(j)=1/(1+ax1);
    ny2(NX-j)=1/(1+ax1);
    ny3(j)=(1-ax1)/(1+ax1);
    ny3(NX-j)=(1-ax1)/(1+ax1);
end

for k=1:pml_length
    ax=(pml_length-k+1)/pml_length;
    ax1=p0*ax^3;
    nz1(k)=ax1;
    nz1(NX-k+1)=ax1;
    mz2(k)=1/(1+ax1);
    mz2(NX-k+1)=1/(1+ax1);
    mz3(k)=(1-ax1)/(1+ax1);
    mz3(NX-k+1)=(1-ax1)/(1+ax1);
    ax=(pml_length-k+0.5)/pml_length;
    ax1=p0*ax^3;
    mz1(k)=ax1;
    mz1(NX-k)=ax1;
    nz2(k)=1/(1+ax1);
    nz2(NX-k)=1/(1+ax1);
    nz3(k)=(1-ax1)/(1+ax1);
    nz3(NX-k)=(1-ax1)/(1+ax1);
end


Ts=0.0;         %time sample
Tsmax=2048;     %maximum number of time samples
Is=0;

% Initialize the maximum matrix for Hx and Hy
Hx_max = zeros(NX, NY);
Hy_max = zeros(NX, NY);
Hz_max = zeros(NX, NY);
%***************************start of FDTD loop******************
for n=1:Tsmax
    Ts=Ts+1;
%***************************Dx field****************************
    for i=2:px
        for j=2:NY
            for k=2:NZ
                curl_h=(hz(i,j,k)-hz(i,j-1,k)-hy(i,j,k)+hy(i,j,k-1));
                curr_dxl(i,j,k)=curr_dxl(i,j,k)+curl_h;
                dx(i,j,k)=(my3(j)*mz3(k)*dx(i,j,k))+(my2(j)*mz2(k)*0.5*(curl_h+mx1(i)*curr_dxl(i,j,k)));
            end
        end
    end
 %***************************************************************
 
    for i=px+1:pxb+1
        for j=2:NY
            for k=2:NZ
                curl_h=(hz(i,j,k)-hz(i,j-1,k)-hy(i,j,k)+hy(i,j,k-1));
                dx(i,j,k)=(my3(j)*mz3(k)*dx(i,j,k))+(my2(j)*mz2(k)*0.5*curl_h);
            end
        end
    end
 %***************************************************************
 for i=pxb+2:NX
     for j=2:NY
         for k=2:NZ
             x_up=i-pxb-1;
             curl_h=(hz(i,j,k)-hz(i,j-1,k)-hy(i,j,k)+hy(i,j,k-1));
             curr_dxh(x_up,j,k)=curr_dxh(x_up,j,k)+curl_h;
             dx(i,j,k)=(my3(j)*mz3(k)*dx(i,j,k))+(my2(j)*mz2(k)*0.5*(curl_h+mx1(i)*curr_dxh(x_up,j,k)));
         end
     end
 end
 
 %************************Dy field*******************************
 for i=2:NX
        for j=2:py
            for k=2:NZ
                curl_h=(hx(i,j,k)-hx(i,j,k-1)-hz(i,j,k)+hz(i-1,j,k));
                curr_dyl(i,j,k)=curr_dyl(i,j,k)+curl_h;
                dy(i,j,k)=(mx3(i)*mz3(k)*dy(i,j,k))+(mx2(i)*mz2(k)*0.5*(curl_h+my1(j)*curr_dyl(i,j,k)));
            end
        end
 end
 %***************************************************************
 for i=2:NX
        for j=py+1:pyb+1
            for k=2:NZ
                curl_h=(hx(i,j,k)-hx(i,j,k-1)-hz(i,j,k)+hz(i-1,j,k));
                dy(i,j,k)=(mx3(i)*mz3(k)*dy(i,j,k))+(mx2(i)*mz2(k)*0.5*curl_h);
            end
        end
 end
 %***************************************************************
 for i=2:NX
        for j=pyb+2:NY
            for k=2:NZ
                y_up=j-pyb-1;
                curl_h=(hx(i,j,k)-hx(i,j,k-1)-hz(i,j,k)+hz(i-1,j,k));
                curr_dyh(i,y_up,k)=curr_dyh(i,y_up,k)+curl_h;
                dy(i,j,k)=(mx3(i)*mz3(k)*dy(i,j,k))+(mx2(i)*mz2(k)*0.5*(curl_h+my1(j)*curr_dyh(i,y_up,k)));
            end
        end
 end
 %***************************Dz field****************************
  for i=2:NX
        for j=2:NY
            for k=2:pz
                curl_h=(hy(i,j,k)-hy(i-1,j,k)-hx(i,j,k)+hx(i,j-1,k));
                curr_dzl(i,j,k)=curr_dzl(i,j,k)+curl_h;
                dz(i,j,k)=(mx3(i)*my3(j)*dz(i,j,k))+(mx2(i)*my2(j)*0.5*(curl_h+mz1(k)*curr_dzl(i,j,k)));
            end
        end
 end
 %***************************************************************
 for i=2:NX
        for j=2:NY
            for k=pz+1:pzb+1
                curl_h=(hy(i,j,k)-hy(i-1,j,k)-hx(i,j,k)+hx(i,j-1,k));
                dz(i,j,k)=(mx3(i)*my3(j)*dz(i,j,k))+(mx2(i)*my2(j)*0.5*curl_h);
            end
        end
 end
 %***************************************************************
 for i=2:NX
        for j=2:NY
            for k=pzb+2:NZ
                z_up=k-pzb-1;
                curl_h=(hy(i,j,k)-hy(i-1,j,k)-hx(i,j,k)+hx(i,j-1,k));
                curr_dzh(i,j,z_up)=curr_dzh(i,j,z_up)+curl_h;
                dz(i,j,k)=(mx3(i)*my3(j)*dz(i,j,k))+(mx2(i)*my2(j)*0.5*(curl_h+mz1(k)*curr_dzh(i,j,z_up)));
            end
        end
 end
 %***********************Source**********************************
 t0=20.0;
 tw=6.0;
 source(n)=exp(-0.5*(((Ts-t0)/tw)^2));   
 dz(xc,yc,zc)=source(n);

 %************************Calculating E from D*******************
  
 % Update the electric field inside the sphere region
for k = 2:NZ-1
    for j = 2:NY-1
        for i = 2:NX-1
            % Calculate the distance to the center of the sphere
            distance_to_center = sqrt((i - load_center(1))^2 + (j - load_center(2))^2 + (k - load_center(3))^2);

            if distance_to_center <= load_radius
                % If the point is inside the sphere, update I and electric field in each direction

                % Update I_x and E_x
                ex(i,j,k) = cx(i,j,k) * (dx(i,j,k) - I_x(i,j,k));  % Update E_x
                I_x(i,j,k) = I_x(i,j,k) + ex(i,j,k);  % Accumulate electric field history effect (x direction)

                % Update I_y and E_y
                ey(i,j,k) = cy(i,j,k) * (dy(i,j,k) - I_y(i,j,k));  % Update E_y
                I_y(i,j,k) = I_y(i,j,k) + ey(i,j,k);  % Accumulate electric field history effect (y direction)

                % Update I_z and E_z
                ez(i,j,k) = cz(i,j,k) * (dz(i,j,k) - I_z(i,j,k));  % Update E_z
                I_z(i,j,k) = I_z(i,j,k) + ez(i,j,k);  % Accumulate electric field history effect (z direction)
            else
                % If the point is outside the sphere, use the default electric field update
                ex(i,j,k) = cx(i,j,k) * dx(i,j,k);
                ey(i,j,k) = cy(i,j,k) * dy(i,j,k);
                ez(i,j,k) = cz(i,j,k) * dz(i,j,k);
            end
        end
    end
end
 %***************************************************************
 for j=1:NY
     for k=1:NZ
         ex(1,j,k)=0;
         ey(1,j,k)=0;
         ez(1,j,k)=0;
     end
 end
 for j=1:NY
     for k=1:NZ
         ex(NX,j,k)=0;
         ey(NX,j,k)=0;
         ez(NX,j,k)=0;
     end
 end
 for i=1:NX
     for k=1:NZ
         ex(i,1,k)=0;
         ey(i,1,k)=0;
         ez(i,1,k)=0;
     end
 end
 for i=1:NX
     for k=1:NZ
         ex(i,NY,k)=0;
         ey(i,NY,k)=0;
         ez(i,NY,k)=0;
     end
 end
 for i=1:NX
     for j=1:NY
         ex(i,j,1)=0;
         ey(i,j,1)=0;
         ez(i,j,1)=0;
     end
 end
 for i=1:NX
     for j=1:NY
         ex(i,j,NZ)=0;
         ey(i,j,NZ)=0;
         ez(i,j,NZ)=0;
     end
 end
 %***************************Hx field****************************
    for i=1:px
        for j=1:NY-1
            for k=1:NZ-1
                curl_e=(ey(i,j,k+1)-ey(i,j,k)-ez(i,j+1,k)+ez(i,j,k));
                curr_hxl(i,j,k)=curr_hxl(i,j,k)+curl_e;
                hx(i,j,k)=(ny3(j)*nz3(k)*hx(i,j,k))+ny2(j)*nz2(k)*0.5*(curl_e+nx1(i)*curr_hxl(i,j,k));
            end
        end
    end
 %***************************************************************
 for i=px+1:pxb+1
        for j=1:NY-1
            for k=1:NZ-1
                curl_e=(ey(i,j,k+1)-ey(i,j,k)-ez(i,j+1,k)+ez(i,j,k));
                hx(i,j,k)=(ny3(j)*nz3(k)*hx(i,j,k))+ny2(j)*nz2(k)*0.5*curl_e;
            end
        end
 end
 %***************************************************************
 for i=pxb+2:NX
        for j=1:NY-1
            for k=1:NZ-1
                x_up=i-pxb-1;
                curl_e=(ey(i,j,k+1)-ey(i,j,k)-ez(i,j+1,k)+ez(i,j,k));
                curr_hxh(x_up,j,k)=curr_hxh(x_up,j,k)+curl_e;
                hx(i,j,k)=(ny3(j)*nz3(k)*hx(i,j,k))+ny2(j)*nz2(k)*0.5*(curl_e+nx1(i)*curr_hxh(x_up,j,k));
            end
        end
 end
 
 %***************************Hy field****************************
    for i=1:NX-1
        for j=1:py
            for k=1:NZ-1
                curl_e=(ez(i+1,j,k)-ez(i,j,k)-ex(i,j,k+1)+ex(i,j,k));
                curr_hyl(i,j,k)=curr_hyl(i,j,k)+curl_e;
                hy(i,j,k)=(nx3(i)*nz3(k)*hy(i,j,k))+nx2(i)*nz2(k)*0.5*(curl_e+ny1(j)*curr_hyl(i,j,k));
            end
        end
    end
 %***************************************************************
 for i=1:NX-1
        for j=py+1:pyb+1
            for k=1:NZ-1
                curl_e=(ez(i+1,j,k)-ez(i,j,k)-ex(i,j,k+1)+ex(i,j,k));
                hy(i,j,k)=(nx3(i)*nz3(k)*hy(i,j,k))+nx2(i)*nz2(k)*0.5*curl_e;
            end
        end
 end
 %***************************************************************
 for i=1:NX-1
        for j=pyb+2:NY
            for k=1:NZ-1
                y_up=j-pyb-1;
                curl_e=(ez(i+1,j,k)-ez(i,j,k)-ex(i,j,k+1)+ex(i,j,k));
                curr_hyh(i,y_up,k)=curr_hyh(i,y_up,k)+curl_e;
                hy(i,j,k)=(nx3(i)*nz3(k)*hy(i,j,k))+nx2(i)*nz2(k)*0.5*(curl_e+ny1(j)*curr_hyh(i,y_up,k));
            end
        end
 end
 %***************************Hz field****************************
    for i=1:NX-1
        for j=1:NY-1
            for k=1:pz
                curl_e=(ex(i,j+1,k)-ex(i,j,k)-ey(i+1,j,k)+ey(i,j,k));
                curr_hzl(i,j,k)=curr_hzl(i,j,k)+curl_e;
                hz(i,j,k)=(nx3(i)*ny3(j)*hz(i,j,k))+nx2(i)*ny2(j)*0.5*(curl_e+nz1(k)*curr_hzl(i,j,k));
            end
        end
    end
 %****************************************************************
 for i=1:NX-1
        for j=1:NY-1
            for k=pz+1:pzb+1
                curl_e=(ex(i,j+1,k)-ex(i,j,k)-ey(i+1,j,k)+ey(i,j,k));
                hz(i,j,k)=(nx3(i)*ny3(j)*hz(i,j,k))+nx2(i)*ny2(j)*0.5*curl_e;
            end
        end
 end
 %****************************************************************
 for i=1:NX-1
        for j=1:NY-1
            for k=pzb+2:NZ
                z_up=k-pzb-1;
                curl_e=(ex(i,j+1,k)-ex(i,j,k)-ey(i+1,j,k)+ey(i,j,k));
                curr_hzh(i,j,z_up)=curr_hzh(i,j,z_up)+curl_e;
                hz(i,j,k)=(nx3(i)*ny3(j)*hz(i,j,k))+nx2(i)*ny2(j)*0.5*(curl_e+nz1(k)*curr_hzh(i,j,z_up));
            end
        end
 end
%*****************************************************************

k1=zc-floor(Nd/2);
k2=zc+floor(Nd/2);
for k=k1:k2
    curr_dipole(k-k1+1,n)=(hx(xc,yc,k)-hx(xc,yc-1,k)-hy(xc,yc,k)+hy(xc-1,yc,k))/Z0;
end
timestep=int2str(n);

% Set the time step range
time_range_start = 200;
time_range_end = 400;

% Compute the sine values during the time steps and extract the maximum value
for n = time_range_start:time_range_end
    % For each time step n, extract Hx and Hy at each point
    for ii = 1:NX
        for jj = 1:NY
            % Extract Hx and Hy components at z=zc plane
            Hx_value = hx(ii, jj, zc);  % Hx component
            Hy_value = hy(ii, jj, zc);  % Hy component

            % Update the maximum value (taking absolute value)
            Hx_max(ii,jj) = max(Hx_max(ii,jj), abs(Hx_value));
            Hy_max(ii,jj) = max(Hy_max(ii,jj), abs(Hy_value));
        end
    end
end
end
% Compute the B1 field
B1 = sqrt(Hx_max.^2 + Hy_max.^2);

% Define a mask to exclude regions outside the sphere
mask = false(NX, NY);
for i = 1:NX
    for j = 1:NY
        % Calculate the distance from each point in the XY plane to the center of the sphere
        d_xy = sqrt((i - load_center(1))^2 + (j - load_center(2))^2);
        % If the distance is within the radius of the sphere, set to true
        if d_xy <= load_radius
            mask(i,j) = true;
        end
    end
end

% Normalize the B1 field only inside the sphere, excluding the region outside the sphere
B1(~mask) = NaN;
B1_normalized = B1 / max(B1(:));

% Plot the normalized B1 field
figure;
surf(1:NX, 1:NY, B1_normalized);
view(0, 90);
shading interp;
colorbar;

% Create a nonlinear colormap
cmap_original = jet(256);  % Use jet colormap
nonlinear_idx = linspace(0, 1, 256).^0.65;  % Increase contrast in low-value regions
cmap_custom = interp1(linspace(0, 1, 256), cmap_original, nonlinear_idx);

% Apply the custom colormap
colormap(cmap_custom);

% Set the colorbar limits
set(gca, 'Clim', [0, 1]);
%-------------------end of FDTD loop---------------------------------------

%-------------------FFT----------------------------------------------------
% Record the current at the center point of the dipole
I_input = curr_dipole(ceil(Nd / 2), :);  % Current at the center of the antenna

% Compute the Fourier Transform (FFT)
I_freq = fft(I_input);

% Get frequency components
n = length(I_input);  % Number of samples
f = (0:n-1) * (1 / (n * dt));  % Frequency vector

% Compute the normalized current magnitude
I_magnitude = abs(I_freq / n);

% Find the frequency corresponding to the maximum component
[~, idx] = max(I_magnitude(1:floor(n / 2)));  % Consider only positive frequencies
resonant_frequency = f(idx);  % Resonant frequency

% Display the resonant frequency
disp(['Resonant frequency is ', num2str(resonant_frequency), ' Hz']);
