function im_out = BFilMap(im, map, radius, sigma_space, sigma_range)

% im_out = im;
% for i = 1+radius : size(im, 1)-radius
%     for j = 1+radius : size(im, 2)-radius
%         if map(i, j) == 1
%             m = i-radius:i+radius;
%             n = j-radius:j+radius;
%             [M, N] = meshgrid(m, n);
%             dist = exp(-((M-i)^2+(N-j)^2)/(2*sigma_space^2));
%             color = exp(-sum((repmat(im(i, j, :), [length(m) length(n) 1])-im(m, n, :)).^2, 3)/(2*sigma_range));
%             summation = sum(sum(im(m, n, :) .* repmat((dist .* color), [1 1 3])));
%             weight = sum(sum(dist .* color));
%             im_out(i, j, :) = summation/weight;
%         end
%     end
% end

im_out = im;
summation = zeros(size(im));
weight = zeros(size(im, 1), size(im, 2));
i = 1:size(im, 1);
j = 1:size(im, 2);
for m = -radius : radius
    for n = -radius : radius
        ii = max(1, min(size(im, 1), i+m));
        jj = max(1, min(size(im, 2), j+n));
        dist = exp(-(m^2+n^2)/(2*sigma_space^2));
        color = exp(-sum((im(i, j, :)-im(ii, jj, :)).^2, 3)/(2*sigma_range));
        summation = summation + repmat(map, [1 1 3]) .* (im(ii, jj, :) .* repmat(dist * color, [1 1 3]));
        weight = weight + map .* (dist .* color);
    end
end
avg = summation ./ repmat(weight, [1 1 3]);
im_out(find(repmat(map, [1 1 3]))) = avg(find(repmat(map, [1 1 3])));

end