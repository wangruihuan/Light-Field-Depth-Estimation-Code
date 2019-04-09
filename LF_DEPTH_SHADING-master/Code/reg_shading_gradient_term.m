function [A_sh, b_sh] = reg_shading_gradient_term(shading, lighting, shading_confi, LF_parameters)

%Parameters
x_size = LF_parameters.x_size;
y_size = LF_parameters.y_size;
num_pixels = x_size * y_size;

confi_diag = spdiags(shading_confi(:),0,num_pixels,num_pixels);

[Gx_sh, Gy_sh] = imgradientxy(shading);

dnx_dx_mat = CONVKERNEL2MATRIX([0 0 0; -2 4 -2; 0 0 0], y_size, x_size);
dny_dx_mat = CONVKERNEL2MATRIX([2 0 -2; 0 0 0; -2 0 2], y_size, x_size);
dnx_dy_mat = CONVKERNEL2MATRIX([2 0 -2; 0 0 0; -2 0 2], y_size, x_size);
dny_dy_mat = CONVKERNEL2MATRIX([0 0 0; -2 4 -2; 0 0 0]', y_size, x_size);

zero_mat = sparse([], [], [], num_pixels, num_pixels);
lx_mat = spdiags(lighting(1)*ones(num_pixels,1), 0, num_pixels, num_pixels);
ly_mat = spdiags(lighting(2)*ones(num_pixels,1), 0, num_pixels, num_pixels);

l_dot_mat = [lx_mat ly_mat zero_mat zero_mat; zero_mat zero_mat lx_mat ly_mat];
partials_mat = [dnx_dx_mat; dny_dx_mat; dnx_dy_mat; dny_dy_mat];

%Outputs
A_sh = [confi_diag zero_mat; zero_mat confi_diag] * l_dot_mat * partials_mat;
b_sh = [confi_diag zero_mat; zero_mat confi_diag] * [Gx_sh(:); Gy_sh(:)];

end