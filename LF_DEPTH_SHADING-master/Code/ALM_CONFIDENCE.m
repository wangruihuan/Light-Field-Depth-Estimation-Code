function confidence = ALM_CONFIDENCE(response,IM_Pinhole,LF_parameters,minmax)

x_size = LF_parameters.x_size;
y_size = LF_parameters.y_size;

luminance = rgb2ycbcr(IM_Pinhole);
luminance = luminance(:,:,1);

std_alm = 0.2;

if(minmax == 1)
    %Defocus
    for x = 1:x_size
        for y = 1:y_size
            response_squeeze = squeeze(response(y,x,:));
            C_1 = max(response_squeeze);
            confidence(y,x) = (1/(sum(exp(-((response_squeeze-C_1).^2/(2*(std_alm^2))))))) * (max(response_squeeze)-min(response_squeeze));
        end
    end
    
else
    %Correspondence
    for x = 1:x_size
        for y = 1:y_size
            response_squeeze = squeeze(response(y,x,:));
            C_1 = min(response_squeeze);
            confidence(y,x) = (1/(sum(exp(-((response_squeeze-C_1).^2/(2*(std_alm^2))))))) * (max(response_squeeze)-min(response_squeeze));
        end
    end
end

end
