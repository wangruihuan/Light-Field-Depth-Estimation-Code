function [A_d, b_d] = reg_data_term_combined(combined_depth, combined_confi)

combined_confi(~isfinite(combined_confi)) = 0;

%Confidence matrix
confi_diag = spdiags(combined_confi(:),0,length(combined_confi(:)),length(combined_confi(:)));

%Outputs
A_d = confi_diag;
b_d = combined_confi(:).*combined_depth(:);

end