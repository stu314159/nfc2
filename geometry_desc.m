clear
clc
close 'all'

%% overall channel dimensions
Lx_p = 1;
Ly_p = 1;
Lz_p = 8;

Ny_divs = 15;

plot_geom = true;

%% select fluid for the simulation
fluid = 1;
% 1 = glycerin
% 2 = glycol
% 3 = water

switch fluid
    case 1
        rho_p = 1260;
        nu_p = 1.49/rho_p;
        
    case 2
        rho_p = 965.3;
        nu_p = 0.06/rho_p;
        
    case 3
        rho_p = 1000;
        nu_p = 1e-3/rho_p;
        
end

obstacle = 3;
% 0 = no obstacle
% 1 = cylinder, bottom to top
% 2 = cylinder, partial height
% 3 = cylindrical piling with elliptical scour pit

switch obstacle
    case 0
        Lo = Ly_p;
        
    case 1
       x_c = 0.5*Lx_p;
       z_c = 0.5*Lz_p;
       cyl_rad = 0.1*Lx_p;
       Lo = cyl_rad*2;
       
    case 2
       x_c = 0.5*Lx_p;
       z_c = 0.5*Lz_p;
       y_max = 0.75*Ly_p;
       cyl_rad = 0.1*Lx_p;
       Lo = cyl_rad*2;

    case 3
       x_c = 0.5*Lx_p; %Center of channel (width)
       z_c = 0.5*Lz_p; %Center of channel (length)
       cyl_rad = 0.1*Lx_p; %Piling radius
       Lo = cyl_rad*2;

	%Scour Pit Dimensions
	ellip_a = 2*2*cyl_rad; %2x diameter of pile
	ellip_b = 1*2*cyl_rad; %1x diam of pile
	ellip_c = 4*2*cyl_rad; %4x diam of pile
	ellip_x = x_c; %center of pile (width)
	ellip_z = z_c + cyl_rad; %center located at back of pile (length)
	ellip_y = ellip_b; %bottom of ellipse at ym
          
    case 4
       x_c = 0.5*Lx_p;
       y_c = 0.5*Ly_p;
       z_c = 0.75*Lz_p;
       sphere_rad = 0.1*Ly_p;
       Lo = sphere_rad*2;
      
end

%% generate the lattice discretization
xm = 0; xp = Lx_p;
ym = 0; yp = Ly_p;
zm = 0; zp = Lz_p;

Ny = ceil((Ny_divs-1)*(Ly_p/Lo))+1;
Nx = ceil((Ny_divs-1)*(Lx_p/Lo))+1;
Nz = ceil((Ny_divs-1)*(Lz_p/Lo))+1;

[gcoord,~,faces]=Brick3Dr2(xm,xp,ym,yp,zm,zp,Nx,Ny,Nz);
[nnodes,~]=size(gcoord);

%% get wall solid nodes
snl = [faces.zx_m; faces.zx_p; faces.zy_m; faces.zy_p]; snl = snl(:);
snl = unique(snl);

%% get inlet nodes
inl = faces.xy_m; 
inl = setxor(inl,intersect(inl,snl)); % eliminate solid nodes from inl

%% get outlet nodes
onl = faces.xy_p;
onl = setxor(onl,intersect(onl,snl)); %eliminate solid nodes from onl

%% find obstacle nodes
switch obstacle
    case 0
        obst_list = [ ];
    case 1
       
       obst_list = find(((gcoord(:,1) - x_c).^2 + (gcoord(:,3)-z_c).^2) < cyl_rad*cyl_rad);

    case 2
      
       obst_list = find((((gcoord(:,1) - x_c).^2 + (gcoord(:,3)-z_c).^2) < cyl_rad*cyl_rad) & ...
         (gcoord(:,2) < y_max));

    case 3

       %% get solid nodes...bottom is thick, top is single plate
       snl = find(gcoord(:,2)<ellip_b); snl =snl(:);
       snl = unique(snl);

       %% get inlet nodes
       inl = faces.xy_m; 
       inl = setxor(inl,intersect(inl,snl)); % eliminate solid nodes from inl

       %% get outlet nodes
	onl = faces.xy_p;
	onl = setxor(onl,intersect(onl,snl)); %eliminate solid nodes from onl

       scourpit = find(((gcoord(:,1) - ellip_x).^2/ellip_a.^2) + ((gcoord(:,2)-ellip_y).^2/ellip_b.^2) + ((gcoord(:,3) - ellip_z).^2/ellip_c.^2) <= 1);

      %Remove ellipsoid nodes from bottom of channel.
      snl = setxor(snl,intersect(snl,scourpit));

      % add the cylindrical obstacle to the obstacle list
      obst_list = find(((gcoord(:,1) - x_c).^2 + (gcoord(:,3)-z_c).^2) < cyl_rad*cyl_rad);
     
    case 4
        
        obst_list = find(((gcoord(:,1) - x_c).^2 + (gcoord(:,2) - y_c).^2 + ...
            (gcoord(:,3) - z_c).^2 < sphere_rad*sphere_rad));
end

% add obstacle nodes to the solid node list
snl = union(snl,obst_list);


%% find pressure reference lattice point
dx = 1/(Ny_divs-1);
l_conv_fact = (dx*Lo);
eps_l = l_conv_fact;
x_pref = 0.5*Lx_p;
y_pref = 0.5*Ly_p;
z_pref = 0.96*Lz_p;
p_ref_LP=find((abs(gcoord(:,1)-x_pref)<=(eps_l/2)) & (abs(gcoord(:,2)-y_pref)<=(eps_l/2)) & ...
    (abs(gcoord(:,3)-z_pref)<=(eps_l/2)));
if(~isempty(p_ref_LP))
    p_ref_LP=p_ref_LP(1);
else
    error('No Reference Pressure Point!!');
end

%% plot the relevant lattice points to confirm correctness
if plot_geom
    figure(1)
%     scatter3(gcoord(inl,1),gcoord(inl,2),gcoord(inl,3),'r.');
%     hold on
%     scatter3(gcoord(onl,1),gcoord(onl,2),gcoord(onl,3),'b.');
     scatter3(gcoord(snl,1),gcoord(snl,2),gcoord(snl,3),'g.'); hold on;
%     hold off
    scatter3(gcoord(obst_list,1),gcoord(obst_list,2),gcoord(obst_list,3),'b.');
     axis([0 Lx_p 0 Ly_p 0 Lz_p]);
    axis equal
    view([-99 52]);
   
end

%% save the data to a *.mat file
file_name = 'geometry_description.mat';
save(file_name,'Lx_p','Ly_p','Lz_p','Lo','Ny_divs','rho_p','nu_p','snl','inl','onl','p_ref_LP');


