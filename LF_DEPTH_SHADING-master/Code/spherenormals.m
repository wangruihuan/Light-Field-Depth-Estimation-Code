%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;
num_pixels = y_size * x_size;

sphere_normals = zeros(y_size, x_size, 3);
[X,Y] = meshgrid(1:x_size,1:y_size);

radius = 170;

for i = 1:x_size
    for j = 1:y_size
        x_ind = X(j,i)-(x_size/2);
        y_ind = Y(j,i)-(y_size/2);
        z_ind = sqrt((radius^2)-(x_ind^2)-(y_ind^2));
        dist = sqrt((x_ind^2)+(y_ind^2));
        if (dist<=radius)
            sphere_normals(j,i,1) = x_ind/radius;
            sphere_normals(j,i,2) = y_ind/radius;
            sphere_normals(j,i,3) = z_ind/radius;
        end
    end
end

