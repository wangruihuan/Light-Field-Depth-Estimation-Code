function combined_response = COMBINE_RESPONSES(defocus_response, corresp_response, defocus_confidence, corresp_confidence)

radius = 1;

combined_response = zeros(size(defocus_response));

for y = 1:size(defocus_response,1)
    for x = 1:size(defocus_response,2)
        weight_sum = 0;
        for j = -radius:radius
            for i = -radius:radius
                y_cord = y+j;
                x_cord = x+i;
                if (y_cord < 1)
                    y_cord = 1;
                elseif (y_cord > size(defocus_response,1))
                    y_cord = size(defocus_response,1);
                end
                if (x_cord < 1)
                    x_cord = 1;
                elseif (x_cord > size(defocus_response,2))
                    x_cord = size(defocus_response,2);
                end
                defocus_response_squeeze = squeeze(defocus_response(y_cord,x_cord,:));
                corresp_response_squeeze = squeeze(corresp_response(y_cord,x_cord,:));
                combined_response(y,x,:) = squeeze(combined_response(y,x,:)) + ...
                    (defocus_confidence(y_cord, x_cord) * defocus_response_squeeze) + ...
                    (corresp_confidence(y_cord, x_cord) * corresp_response_squeeze);
                weight_sum = weight_sum + defocus_confidence(y_cord, x_cord) + corresp_confidence(y_cord, x_cord);
            end
        end
        combined_response(y,x,:) = combined_response(y,x,:)/weight_sum;
    end
end

end