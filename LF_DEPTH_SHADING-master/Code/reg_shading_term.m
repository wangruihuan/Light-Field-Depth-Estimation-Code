function [A_sh, b_sh] = reg_shading_term(shading, lighting, shading_confi, LF_parameters)

%Parameters
x_size = LF_parameters.x_size;
y_size = LF_parameters.y_size;
num_pixels = x_size * y_size;

confi_diag = spdiags(shading_confi(:),0,num_pixels,num_pixels);

nx_mat = CONVKERNEL2MATRIX([0 0 0; 2 0 -2; 0 0 0], y_size, x_size);
ny_mat = CONVKERNEL2MATRIX([0 -2 0; 0 0 0; 0 2 0], y_size, x_size);

lx_mat = spdiags(lighting(1)*ones(num_pixels,1), 0, num_pixels, num_pixels);
ly_mat = spdiags(lighting(2)*ones(num_pixels,1), 0, num_pixels, num_pixels);

l_dot_mat = [lx_mat ly_mat];
normals_mat = [nx_mat; ny_mat];

%Outputs
A_sh = confi_diag * l_dot_mat * normals_mat;
b_sh = confi_diag * (shading(:) - (lighting(3) * 4));

end