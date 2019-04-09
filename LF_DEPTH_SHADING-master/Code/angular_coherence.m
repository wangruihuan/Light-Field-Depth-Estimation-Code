function [Angular_Coherence_Mat, angular_use_indices] = angular_coherence(initial_depth, LF_Remap, LF_parameters)

LF_y_size = LF_parameters.LF_y_size;
LF_x_size = LF_parameters.LF_x_size;
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;
UV_diameter = LF_parameters.UV_diameter;

LF_Focused_Remap = zeros(LF_y_size,LF_x_size,3);
LF_Remap_Indices = zeros(LF_y_size,LF_x_size,8);
Angular_Coherence_ii = zeros(UV_diameter*UV_diameter*4*x_size*y_size,1);
Angular_Coherence_jj = zeros(UV_diameter*UV_diameter*4*x_size*y_size,1);
Angular_Coherence_ss = zeros(UV_diameter*UV_diameter*4*x_size*y_size,1);
count = 1;
for y = 1:y_size
    for x = 1:x_size
        depth   = initial_depth(y,x);
        angular  = get_angular(x,y,depth,LF_Remap, LF_parameters);
        indices = get_inverse_angular_indices(x,y,depth, LF_parameters);
        
        x_min   = (x-1)*UV_diameter+1;
        x_max   = (x)  *UV_diameter;
        y_min   = (y-1)*UV_diameter+1;
        y_max   = (y)  *UV_diameter;
        
        LF_Focused_Remap(y_min:y_max,x_min:x_max,1) = angular(:,:,1);
        LF_Focused_Remap(y_min:y_max,x_min:x_max,2) = angular(:,:,2);
        LF_Focused_Remap(y_min:y_max,x_min:x_max,3) = angular(:,:,3);
        
        LF_Remap_Indices(y_min:y_max,x_min:x_max,1) = indices(:,:,1);
        LF_Remap_Indices(y_min:y_max,x_min:x_max,2) = indices(:,:,2);
        LF_Remap_Indices(y_min:y_max,x_min:x_max,3) = indices(:,:,3);
        LF_Remap_Indices(y_min:y_max,x_min:x_max,4) = indices(:,:,4);
        LF_Remap_Indices(y_min:y_max,x_min:x_max,5) = indices(:,:,5);
        LF_Remap_Indices(y_min:y_max,x_min:x_max,6) = indices(:,:,6);
        LF_Remap_Indices(y_min:y_max,x_min:x_max,7) = indices(:,:,7);
        LF_Remap_Indices(y_min:y_max,x_min:x_max,8) = indices(:,:,8);
        
        spatial_index = ((x-1)*y_size) + y;
        for s = 1:UV_diameter
            for t = 1:UV_diameter
                angular_index = ((t-1)*UV_diameter)+s;
                angular_offset = ((angular_index-1)*x_size*y_size);
                x_1 = indices(s,t,1);
                y_1 = indices(s,t,2);
                x_2 = indices(s,t,3);
                y_2 = indices(s,t,4);
                x_1_spatial = ceil(x_1/UV_diameter);
                y_1_spatial = ceil(y_1/UV_diameter);
                x_2_spatial = ceil(x_2/UV_diameter);
                y_2_spatial = ceil(y_2/UV_diameter);
                pt_1 = ((x_1_spatial-1)*y_size) + y_1_spatial;
                pt_2 = ((x_1_spatial-1)*y_size) + y_2_spatial;
                pt_3 = ((x_2_spatial-1)*y_size) + y_1_spatial;
                pt_4 = ((x_2_spatial-1)*y_size) + y_2_spatial;
                Angular_Coherence_ii(count:count+3) = spatial_index;
                Angular_Coherence_jj(count:count+3) = [pt_1 pt_2 pt_3 pt_4]+angular_offset;
                if (s==4 && t==4)
                    Angular_Coherence_ss(count:count+3) = 1;
                else 
                    Angular_Coherence_ss(count:count+3) = -1;
                end
                count = count+4;
            end
        end
    end
end

Angular_Coherence_Mat = sparse(Angular_Coherence_ii, Angular_Coherence_jj, Angular_Coherence_ss, x_size*y_size, x_size*y_size*UV_diameter*UV_diameter);

angular_use_indices = [24 25];
angular_use_offsets = ((angular_use_indices-1)*x_size*y_size);
keep_indices = zeros(1,x_size*y_size*UV_diameter*UV_diameter);
for r = 1:length(angular_use_indices)
    keep_indices(angular_use_offsets(r)+1:angular_use_offsets(r)+(x_size*y_size)) = r;
end
Angular_Coherence_Mat(:,keep_indices==0) = [];