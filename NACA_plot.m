function [upperSurface,lowerSurface] = NACA_plot(foilDescriptor,z_values,chordlength)
%NACA_plot takes the descriptor of any four-digit NACA foil and creates
%a coordinate file.  
%   INPUTS: foilDescriptor

disp('Initiating NACA_plot')

% first, check to make sure foil_descriptor is a string
if ischar(foilDescriptor) == 0
    disp('This function requires you to enter the input argument')
    disp('as a string variable (i.e. in single quotes).')
    disp('Please try again.')
    return %exits the p
end
% next, pull the required information from the foil descriptor
m = str2double(foilDescriptor(1))/100;
p = str2double(foilDescriptor(2))/10;
t = str2double(foilDescriptor(3:end))/100;

% create the x-vector
x = z_values' - min(z_values); %create vector of x values starting from 0
c = chordlength;
n = length(z_values);

% calculate the y-coordinates for a symmetric airfoil of thickness t
y_t = t*c/0.2* ...
    ((0.2969)*sqrt(x/c)+ ...
    (-0.1260)*(x/c).^1+ ...
    (-0.3516)*(x/c).^2+ ...
    ( 0.2843)*(x/c).^3+ ...
    (-0.1015)*(x/c).^4);

% distinguish between cambered and non-cambered airfoils
if m == 0 % If the airfoil is symmetric, 
    upperSurface = [x' y_t'];
    lowerSurface = [x' -1*y_t'];
else  %if the airfoil is cambered
    
    % use a for loop with conditional logic to calculate the camber line
    for k = 1:length(x)
        if x(k) <= p*c
            y_c(k) = m*(x(k)/p^2).*(2*p-(x(k)/c)); 
            dycdx(k) = 2*m/p^2*(p-(x(k)/c));
        else
            y_c(k) = m*((c-x(k))/(1-p)^2).*(1+(x(k)/c)-2*p);
            dycdx(k) = 2*m/(1-p)^2*(p-(x(k)/c));
        end
    end
    
    % Calculate the coordinates for the upper and lower surfaces
    theta = atan(dycdx);
    x_u = x   - y_t.*sin(theta);
    x_l = x   + y_t.*sin(theta);
    y_u = y_c + y_t.*cos(theta);
    y_l = y_c - y_t.*cos(theta);
    
    % The way the equations work out, sometimes the chord length extends a
    % little beyond what it should.  Rescale so that it's exactly what was
    % entered above.
    scalingFactor = c/(max(x_u)-min(x_u));
    x_u = x_u.*scalingFactor;
    x_l = x_l.*scalingFactor;
    y_u = y_u.*scalingFactor;
    y_l = y_l.*scalingFactor;
    
    upperSurface = [x_u' y_u'];
    lowerSurface = [x_l' y_l'];
        
end


% ensure that the airfoil is closed by setting the leading edge to 
% 0,0 and the trailing edge is at c,0
upperSurface(1,:) = [0 0];
upperSurface(end,:) = [c 0];
lowerSurface(1,:) = [0 0];
lowerSurface(end,:) = [c 0];

% Create a z-vector of zeros (required by SolidWorks) and append to the
% airfoil coordinate files
%upperSurface = [upperSurface zeros(length(upperSurface(:,1)),1)];
%lowerSurface = [lowerSurface zeros(length(lowerSurface(:,1)),1)];



end

