function im_out = RemapAverage(im, sigma, LF_parameters)

UV_diameter = LF_parameters.UV_diameter;
i = 1:LF_parameters.y_size;
j = 1:LF_parameters.x_size;
im_out = zeros(i(end), j(end), 3);
weight = 0;
mid = ceil(UV_diameter/2);
for m = 1:UV_diameter
    for n = 1:UV_diameter
         im_out(i, j, :) = im_out(i, j, :) + exp(-((m-mid)^2+(n-mid)^2)/(2*sigma^2)) ...
             * im((i-1)*UV_diameter+m, (j-1)*UV_diameter+n, :);
         weight = weight + exp(-((m-mid)^2+(n-mid)^2)/(2*sigma^2));
    end
end
im_out = im_out / weight;

end