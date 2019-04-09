function [A_sh, b_sh] = reg_shading_projection_term(normals, LF_parameters)

%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;
num_pixels = y_size * x_size;

%Matrix to reproject pixels into scene
depth_to_points_matrix = DEPTH2POINTS_MATRIX(LF_parameters);

%Matrix to perform dot product for shading constraint
dot_1_diagx = zeros(num_pixels,1);
dot_1_diagy = zeros(num_pixels,1);
dot_1_diagz = zeros(num_pixels,1);
dot_2_diagx = zeros(num_pixels,1);
dot_2_diagy = zeros(num_pixels,1);
dot_2_diagz = zeros(num_pixels,1);
dot_3_diagx = zeros(num_pixels,1);
dot_3_diagy = zeros(num_pixels,1);
dot_3_diagz = zeros(num_pixels,1);
dot_4_diagx = zeros(num_pixels,1);
dot_4_diagy = zeros(num_pixels,1);
dot_4_diagz = zeros(num_pixels,1);

for i = 1:x_size %Calculate sum of normals and modify dot product matrices
    for j = 1:y_size
        i_cord1 = i - 1;
        i_cord2 = i + 1;
        j_cord3 = j - 1;
        j_cord4 = j + 1;
        if (i==1)
            i_cord1 = i;
        elseif (i==x_size)
            i_cord2 = i;
        end
        if (j==1)
            j_cord3 = j;
        elseif (j==y_size)
            j_cord4 = j;
        end
        nsum1 = squeeze(normals(j,i,:) + normals(j,i_cord1,:)); %Normal above
        nsum2 = squeeze(normals(j,i,:) + normals(j,i_cord2,:)); %Normal below
        nsum3 = squeeze(normals(j,i,:) + normals(j_cord3,i,:)); %Normal to left
        nsum4 = squeeze(normals(j,i,:) + normals(j_cord4,i,:)); %Normal to right
        nsum1 = nsum1 / norm(nsum1);
        nsum2 = nsum2 / norm(nsum2);
        nsum3 = nsum3 / norm(nsum3);
        nsum4 = nsum4 / norm(nsum4);
        
        index = j + ((i-1) * y_size);
        dot_1_diagx(index,1) = nsum1(1,1);
        dot_1_diagy(index,1) = nsum1(2,1);
        dot_1_diagz(index,1) = nsum1(3,1);
        dot_2_diagx(index,1) = nsum2(1,1);
        dot_2_diagy(index,1) = nsum2(2,1);
        dot_2_diagz(index,1) = nsum2(3,1);
        dot_3_diagx(index,1) = nsum3(1,1);
        dot_3_diagy(index,1) = nsum3(2,1);
        dot_3_diagz(index,1) = nsum3(3,1);
        dot_4_diagx(index,1) = nsum4(1,1);
        dot_4_diagy(index,1) = nsum4(2,1);
        dot_4_diagz(index,1) = nsum4(3,1);
    end
end

diag_vals_1x = [ones(num_pixels,1).*dot_1_diagx -1*ones(num_pixels,1).*dot_1_diagx];
diag_vals_1y = [ones(num_pixels,1).*dot_1_diagy -1*ones(num_pixels,1).*dot_1_diagy];
diag_vals_1z = [ones(num_pixels,1).*dot_1_diagz -1*ones(num_pixels,1).*dot_1_diagz];
diag_vals_2x = [ones(num_pixels,1).*dot_2_diagx -1*ones(num_pixels,1).*dot_2_diagx];
diag_vals_2y = [ones(num_pixels,1).*dot_2_diagy -1*ones(num_pixels,1).*dot_2_diagy];
diag_vals_2z = [ones(num_pixels,1).*dot_2_diagz -1*ones(num_pixels,1).*dot_2_diagz];
diag_vals_3x = [ones(num_pixels,1).*dot_3_diagx -1*ones(num_pixels,1).*dot_3_diagx];
diag_vals_3y = [ones(num_pixels,1).*dot_3_diagy -1*ones(num_pixels,1).*dot_3_diagy];
diag_vals_3z = [ones(num_pixels,1).*dot_3_diagz -1*ones(num_pixels,1).*dot_3_diagz];
diag_vals_4x = [ones(num_pixels,1).*dot_4_diagx -1*ones(num_pixels,1).*dot_4_diagx];
diag_vals_4y = [ones(num_pixels,1).*dot_4_diagy -1*ones(num_pixels,1).*dot_4_diagy];
diag_vals_4z = [ones(num_pixels,1).*dot_4_diagz -1*ones(num_pixels,1).*dot_4_diagz];
dot_1x = spdiags(diag_vals_1x, [0 -1], num_pixels, num_pixels); %Pixel above
dot_1y = spdiags(diag_vals_1y, [0 -1], num_pixels, num_pixels); 
dot_1z = spdiags(diag_vals_1z, [0 -1], num_pixels, num_pixels); 
dot_2x = spdiags(diag_vals_2x, [0 1], num_pixels, num_pixels); %Pixel below
dot_2y = spdiags(diag_vals_2y, [0 1], num_pixels, num_pixels); %
dot_2z = spdiags(diag_vals_2z, [0 1], num_pixels, num_pixels); 
dot_3x = spdiags(diag_vals_3x, [0 -y_size], num_pixels, num_pixels); %Pixel to left
dot_3y = spdiags(diag_vals_3y, [0 -y_size], num_pixels, num_pixels); 
dot_3z = spdiags(diag_vals_3z, [0 -y_size], num_pixels, num_pixels); 
dot_4x = spdiags(diag_vals_4x, [0 y_size], num_pixels, num_pixels); %Pixel to right
dot_4y = spdiags(diag_vals_4y, [0 y_size], num_pixels, num_pixels); 
dot_4z = spdiags(diag_vals_4z, [0 y_size], num_pixels, num_pixels); 

dot_prod_matrix = [dot_1x, dot_1y, dot_1z;...
                   dot_2x, dot_2y, dot_2z;...
                   dot_3x, dot_3y, dot_3z;...
                   dot_4x, dot_4y, dot_4z];
               
A_sh = dot_prod_matrix * depth_to_points_matrix;
b_sh = zeros(num_pixels,1);

end