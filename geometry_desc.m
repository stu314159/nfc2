clear
clc
close 'all'

%% overall channel dimensions
Lx_p = 1;
Ly_p = 2;
Lz_p = 10;

Ny_divs = input('Ny_divs? Must be ODD-numbered \n\n');clc;

plot_geom = true;

%% select fluid for the simulation
fluid = menu('Fluid?','glycerin','glycol','water');
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

obstacle = menu('Which Obstacle?','None','Cylinder from Top to Bottom','Cylinder:Partial Height'...
    ,'Sphere','Cone','Cone With Angle of Attack','NACA Airfoil') - 1;
% 0 = no obstacle
% 1 = cylinder, bottom to top
% 2 = cylinder, partial height
% 3 = sphere
% 4 = cone
% 5 = cone with angle of attack
% 6 = NACA Airfoil

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
       x_c = 0.5*Lx_p;
       y_c = 0.5*Ly_p;
       z_c = 0.75*Lz_p;
       sphere_rad = 0.1*Ly_p;    
       Lo = sphere_rad*2;
       
    case 4
        x_c = 0.5*Lx_p;
        y_c = 0.5*Ly_p;
        r_max = 0.25*Lx_p; r_min = 0; 
        z_min = 0.4*Lz_p;
        z_max = 0.6*Lz_p;
        cone_len = z_max - z_min;
        Lo = cone_len;
        
    case 5
        angle_of_attack = 0;
        x_c = 0.5*Lx_p;
        r_max = 0.25*Lx_p; r_min = 0;
        z_c_max = 0.6*Lz_p-(0.6 - 0.5)*Lz_p*cosd(angle_of_attack);
        z_c_min = 0.5*Lz_p - (0.5*Lz_p - 0.4*Lz_p)*cosd(angle_of_attack);
        cone_len = z_c_max - z_c_min;
        Lo = cone_len;
        y_c_zmin = 0.5*Ly_p + (0.5-0.4)*Lz_p*sind(angle_of_attack);
        y_c_zmax = 0.5*Ly_p - (0.5-0.4)*Lz_p*sind(angle_of_attack);
        
    case 6 
        NACA_foil_designator = input('Type in the FOUR-digit NACA Airfoild Designator \n\n','s');
        z_min = 2;
        z_max = 5;
        x_c = 0.5*Lx_p;
        chordlength = z_max - z_min; % shouldn't be longer than 3 because it has to fit in channel of length 8
        Lo = chordlength;
        
             
end

%% generate the lattice discretization
xm = 0; xp = Lx_p;
ym = 0; yp = Ly_p;
zm = 0; zp = Lz_p;

Ny = ceil((Ny_divs-1)*(Ly_p./Lo))+1;
Nx = ceil((Ny_divs-1)*(Lx_p./Lo))+1;
Nz = ceil((Ny_divs-1)*(Lz_p./Lo))+1;

[gcoord,~,faces]=Brick3Dr2(xm,xp,ym,yp,zm,zp,Nx,Ny,Nz);
[nnodes,~]=size(gcoord);

%% get wall solid nodes
snl = [faces.zx_m; faces.zx_p; faces.zy_m; faces.zy_p]; snl = snl(:);
snl = unique(snl);

%% get in
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
        
        obst_list = find(((gcoord(:,1) - x_c).^2 + (gcoord(:,2) - y_c).^2 + ...
            (gcoord(:,3) - z_c).^2 < sphere_rad*sphere_rad));
        
    case 4 
        
        % z >= z_min; z<= z_max and distance from x_c, y_c within r(z)
        r_z = @(z) (z - z_min).*(r_max/cone_len);
        obst_list = find((gcoord(:,3)>=z_min) & (gcoord(:,3)<=z_max) & ... % z>=z_min and z<=z_max
            ((gcoord(:,1)-x_c).^2 + (gcoord(:,2)-y_c).^2 <= r_z(gcoord(:,3)).^2)); 
        
    case 5
        r_z = @(z) (z - z_c_min).*(r_max/cone_len);%?????
        y_c = @(z) (z - z_c_min).*(y_c_zmax - y_c_zmin)/cone_len + y_c_zmin; %?????
        z_max = @(y) (y).*(tand(angle_of_attack)) + z_c_max - 0.5.*Ly_p.*tand(angle_of_attack);
        
        obst_list = find((gcoord(:,3)>=z_c_min) & (gcoord(:,3)<=z_max(gcoord(:,2))) & ...
            ((gcoord(:,1) - x_c).^2 + (gcoord(:,2) -y_c(gcoord(:,3))).^2<= r_z(gcoord(:,3)).^2));

    case 6
        z_indices = find(gcoord(:,3) >= z_min & gcoord(:,3) <= z_max);
        z_values = gcoord(z_indices,3);
        foilDesignator = num2str(NACA_foil_designator);
        [upperSurface,lowerSurface] = NACA_plot(foilDesignator,z_values,chordlength);
        
        %values of the upper surface
        z_upper = upperSurface(:,1) + z_min;
        y_upper = upperSurface(:,2) + 0.5*Ly_p;
        upperSurface = [z_upper,y_upper]; 
        
        %values of lower surface
        z_lower = lowerSurface(:,1) + z_min;
        y_lower = lowerSurface(:,2) + 0.5*Ly_p;
        lowerSurface = [z_lower, y_lower];

        %replace y_values of gcoord with y_lower and y_upper at z_indices
        gcoordupper = gcoord(:,:);
        gcoordlower = gcoord(:,:);
        gcoordupper(z_indices,2) = y_upper; 
        gcoordlower(z_indices,2) = y_lower;

        x_min = 0;
        x_max = Lx_p;
        obst_list = find((gcoord(:,2) <= gcoordupper(:,2) & gcoord(:,2) >= gcoordlower(:,2)) &...
            (gcoord(:,3) >= z_min & gcoord(:,3) <= z_max));
        
end
% add obstacle nodes to the solid node listy_c = @(z) (z - z_c_min).*(y_c_min/cone_len);
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
     scatter3(gcoord(inl,1),gcoord(inl,2),gcoord(inl,3),'r.');
     hold on
     scatter3(gcoord(onl,1),gcoord(onl,2),gcoord(onl,3),'b.');
  %   scatter3(gcoord(snl,1),gcoord(snl,2),gcoord(snl,3),'g.'); hold on;
  %   hold off
    scatter3(gcoord(obst_list,1),gcoord(obst_list,2),gcoord(obst_list,3),'b.');
    axis([0 Lx_p 0 Ly_p 0 Lz_p]);
    axis equal
    view([-99 52]);
end

%% save the data to a *.mat file
file_name = 'geometry_description.mat';
save(file_name,'Lx_p','Ly_p','Lz_p','Lo','Ny_divs','rho_p','nu_p','snl','inl','onl','p_ref_LP');


