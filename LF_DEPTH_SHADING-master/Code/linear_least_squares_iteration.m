function [output] = linear_least_squares_iteration(A, b, LF_parameters)

%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;
iter_max = LF_parameters.iter_max;
convergence_ratio = LF_parameters.convergence_ratio;
err_epsilon = LF_parameters.err_epsilon;

output = lsqlin(A, b);

%Iterate
mean_err = zeros(iter_max,1);
for i = 1:iter_max
    sqrerr = sqrt(((A * output) - b).^2);
    err_weights = sqrt(1./(sqrerr + err_epsilon));
    err_weights_sparse = sparse(1:length(err_weights),1:length(err_weights),err_weights);
    mean_err(i,1) = mean(sqrerr(:));
    %fprintf('iteration: %i  mean error: %f \n', i, mean_err(i,1));
    if (i>1 && (abs(mean_err(i-1,1) - mean_err(i,1)) / mean_err(i,1))<convergence_ratio)
        fprintf('Converged \n');
        break;
    end
    A_weighted = err_weights_sparse * A;
    b_weighted = err_weights_sparse * b;
    output = lsqlin(A_weighted, b_weighted);
    %output_image = reshape(output, [y_size x_size]);
    %figure; imagesc(output_image); colormap gray;
end

output = reshape(output, [y_size x_size]);

end