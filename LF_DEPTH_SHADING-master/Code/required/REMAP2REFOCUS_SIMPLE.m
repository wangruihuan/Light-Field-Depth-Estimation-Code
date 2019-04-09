function IM_Remap = REMAP2REFOCUS_SIMPLE(LF_Remap)

UV_diameter   = 7                                   ;
UV_radius     = 3                                   ;
x_size        = size(LF_Remap,2)/UV_diameter        ;
y_size        = size(LF_Remap,1)/UV_diameter        ;

IM_Remap = zeros(y_size,x_size,3);

h = fspecial('average',[UV_diameter UV_diameter]);
LF_Remap = imfilter(LF_Remap,h,'symmetric');

for x = 1:x_size
    for y = 1:y_size
        x_coord_center = 1 + UV_radius + UV_diameter * (x - 1)            ;
        y_coord_center = 1 + UV_radius + UV_diameter * (y - 1)            ;
        
        IM_Remap(y,x,1) = LF_Remap(y_coord_center,x_coord_center,1)     ;
        IM_Remap(y,x,2) = LF_Remap(y_coord_center,x_coord_center,2)     ;
        IM_Remap(y,x,3) = LF_Remap(y_coord_center,x_coord_center,3)     ;
    end
end
end