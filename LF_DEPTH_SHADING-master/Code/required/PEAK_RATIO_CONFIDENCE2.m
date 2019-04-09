function confidence = PEAK_RATIO_CONFIDENCE2(response,LF_parameters,minmax)
%PEAK_RATIO_CONFIDENCE
%           Takes the image refocused to alpha, and outputs response
%           for each pixel. This is a rudimentry contrast-based peak detection.
%           Input : response
%                   minmax        (0 for min- corresp, 1 for max - defocus)
%           Output: confidence

%           EQUATION (7) in paper

x_size           = LF_parameters.x_size                                   ;
y_size           = LF_parameters.y_size                                   ;
depth_resolution = LF_parameters.depth_resolution                         ;

radius           = round(0.1*256)                                         ;

if(minmax == 1)
    % defocus
    for x = 1:x_size
        for y = 1:y_size
            response_squeeze = squeeze(response(y,x,:))              ;
            [C_1,I] = max(response_squeeze)                  ;
            I_min   = max(1,(I-radius))              ;
            I_max   = min(depth_resolution, I+radius);
            [C_2,I] = max([response_squeeze(1:I_min); response_squeeze(I_max:depth_resolution)]);
            confidence(y,x) = C_2/C_1                     ;
            
            if (confidence(y,x) > 1 || isinf(confidence(y,x)))
                confidence(y,x) = 1;
            else
                confidence(y,x) = 1 - confidence(y,x);
            end
        end
    end
    
else
    % correspondence
    for x = 1:x_size
        for y = 1:y_size
            response_squeeze = squeeze(response(y,x,:))              ;
            [C_1,I] = min(response_squeeze)                  ;
            I_min   = max(1,(I-radius))              ;
            I_max   = min(depth_resolution, I+radius);
            [C_2,I] = min([response_squeeze(1:I_min); response_squeeze(I_max:depth_resolution)]);
            confidence(y,x) = C_1/C_2        ;
            
            if (confidence(y,x) > 1 || isinf(confidence(y,x)))
                confidence(y,x) = 1;
            else
                confidence(y,x) = 1 - confidence(y,x);
            end
        end
    end
end

end

