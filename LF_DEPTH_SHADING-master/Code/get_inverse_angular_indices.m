function im_out_remap = get_inverse_angular_indices(x,y,depth, LF_parameters)
%GET_ANGULAR Summary of this function goes here
%   Detailed explanation goes here

UV_radius = LF_parameters.UV_radius;
UV_diameter = LF_parameters.UV_diameter;
depth_resolution = LF_parameters.depth_resolution;
alpha_max = LF_parameters.alpha_max;
alpha_min = LF_parameters.alpha_min;
width     = LF_parameters.x_size;
height    = LF_parameters.y_size;

alpha_step       = (alpha_max-alpha_min)/depth_resolution                 ;
alpha            = alpha_min + (depth-1)*alpha_step                       ;

im_out_remap = zeros(UV_diameter,UV_diameter,8);

for i = -UV_radius+1:UV_radius+1
    for j = -UV_radius+1:UV_radius+1
        
        x_ind   = i*(1-1/alpha) + x;
        y_ind   = j*(1-1/alpha) + y;
        
        x_floor = floor(x_ind);
        y_floor = floor(y_ind);
        
        x_1     = max(1,min(x_floor  ,width ));
        y_1     = max(1,min(y_floor  ,height));
        x_2     = max(1,min(x_floor+1,width ));
        y_2     = max(1,min(y_floor+1,height));
        
        x_1_w   = 1-(x_ind-x_floor)        ;
        x_2_w   = 1-x_1_w                  ;
        y_1_w   = 1-(y_ind-y_floor)        ;
        y_2_w   = 1-y_1_w                  ;
        
        x_1_index = i+UV_radius + (x_1-1)*UV_diameter   ;
        y_1_index = j+UV_radius + (y_1-1)*UV_diameter   ;
        x_2_index = i+UV_radius + (x_2-1)*UV_diameter   ;
        y_2_index = j+UV_radius + (y_2-1)*UV_diameter   ;
        
        x_index_remap = i+UV_radius;
        y_index_remap = j+UV_radius;
        
        im_out_remap(y_index_remap,x_index_remap,1) = x_1_index;
        im_out_remap(y_index_remap,x_index_remap,2) = y_1_index;
        im_out_remap(y_index_remap,x_index_remap,3) = x_2_index;
        im_out_remap(y_index_remap,x_index_remap,4) = y_2_index;
        im_out_remap(y_index_remap,x_index_remap,5) = y_1_w*x_1_w;
        im_out_remap(y_index_remap,x_index_remap,6) = y_2_w*x_1_w;
        im_out_remap(y_index_remap,x_index_remap,7) = y_1_w*x_2_w;
        im_out_remap(y_index_remap,x_index_remap,8) = y_2_w*x_2_w;
    end
end
