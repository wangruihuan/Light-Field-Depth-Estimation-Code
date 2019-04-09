shading = shading_regularization_LF(Angular_Coherence_Mat, angular_use_indices, LF_Remap, normals, normal_clusters_mat, texture_clusters_mat, LF_parameters);

%%Lighting Estimation
lighting = LIGHTINGESTIMATION(normals, shading, LF_parameters);

%%Depth Optimization with Shading
depth_sh = depth_regularization_shading(pinhole_filt, initial_depth, combined_confi, shading, lighting, LF_parameters);
