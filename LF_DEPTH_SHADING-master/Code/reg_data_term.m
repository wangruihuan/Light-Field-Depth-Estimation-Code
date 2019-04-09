function [A_d, b_d] = reg_data_term(defocus_depth,corresp_depth,...
    defocus_confi,corresp_confi)

defocus_confi(~isfinite(defocus_confi)) = 0;
corresp_confi(~isfinite(corresp_confi)) = 0;

%Combine defocus/correspondence depths and confidences
confi_ratio = 1;
depth_buffer        = (defocus_confi > corresp_confi*confi_ratio).*defocus_depth...
                        + (defocus_confi <= corresp_confi*confi_ratio).*corresp_depth;
confi_buffer        = max(defocus_confi, corresp_confi);

%Confidence matrix
confi_diag = spdiags(confi_buffer(:),0,length(confi_buffer(:)),length(confi_buffer(:)));

%Outputs
A_d = confi_diag;
b_d = confi_buffer(:).*depth_buffer(:);

end