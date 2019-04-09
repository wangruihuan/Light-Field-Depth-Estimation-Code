function LF_Remap_Out = reconstruct_original(LF_Remap_Indices,LF_Remap_In, LF_paramters)
%RECONSTRUCT_ORIGINAL(LF_REMAP_INDICES,LF_FOCUSED_REMAP, LF_PARAMTERS) Summary of this function goes here
%   Detailed explanation goes here

LF_y_size = LF_paramters.LF_y_size;
LF_x_size = LF_paramters.LF_x_size;

LF_Remap_Out            = zeros(size(LF_Remap_In));
LF_Remap_Out_Weights    = zeros(size(LF_Remap_In,1),size(LF_Remap_In,2));

for y = 1:LF_y_size
    for x = 1:LF_x_size
        x_1_index = LF_Remap_Indices(y,x,1);
        y_1_index = LF_Remap_Indices(y,x,2);
        x_2_index = LF_Remap_Indices(y,x,3);
        y_2_index = LF_Remap_Indices(y,x,4);
        w_1       = LF_Remap_Indices(y,x,5);
        w_2       = LF_Remap_Indices(y,x,6);
        w_3       = LF_Remap_Indices(y,x,7);
        w_4       = LF_Remap_Indices(y,x,8);
        for c = 1:size(LF_Remap_In,3)
            LF_Remap_Out(y_1_index,x_1_index,c) = LF_Remap_Out(y_1_index,x_1_index,c) + w_1*LF_Remap_In(y,x,c);
            LF_Remap_Out(y_2_index,x_1_index,c) = LF_Remap_Out(y_2_index,x_1_index,c) + w_2*LF_Remap_In(y,x,c);
            LF_Remap_Out(y_1_index,x_2_index,c) = LF_Remap_Out(y_1_index,x_2_index,c) + w_3*LF_Remap_In(y,x,c);
            LF_Remap_Out(y_2_index,x_2_index,c) = LF_Remap_Out(y_2_index,x_2_index,c) + w_4*LF_Remap_In(y,x,c);
        end
        LF_Remap_Out_Weights(y_1_index,x_1_index) = LF_Remap_Out_Weights(y_1_index,x_1_index) + w_1;
        LF_Remap_Out_Weights(y_2_index,x_1_index) = LF_Remap_Out_Weights(y_2_index,x_1_index) + w_2;
        LF_Remap_Out_Weights(y_1_index,x_2_index) = LF_Remap_Out_Weights(y_1_index,x_2_index) + w_3;
        LF_Remap_Out_Weights(y_2_index,x_2_index) = LF_Remap_Out_Weights(y_2_index,x_2_index) + w_4;
    end
end

for c = 1:size(LF_Remap_In,3)
    LF_Remap_Out(:,:,c) = LF_Remap_Out(:,:,c)./LF_Remap_Out_Weights;
end

