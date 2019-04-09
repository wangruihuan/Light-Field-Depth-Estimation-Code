h_start_odd  = 13;
h_start_even = 22;
h_advan  = 8.665;

w_start_odd_o  = 9;
w_start_even_o = 4;
w_advan = 10.004;

x_shift_const  = 0.050;
y_shift_const  =  -0.055;

im_in_h = 3280;
im_in_w = 3280;

    % odd
    y_count = -1;

    image_cords = zeros(400, 400, 2);

    for y_o = h_start_odd : h_advan*2 : im_in_h

        y_count = y_count + 2;
        x_count = 1;

        w_start_odd = w_start_odd_o + y_count*x_shift_const;

        for x = w_start_odd : w_advan : im_in_w

            y = y_o + x_count*y_shift_const;
            % image
            x_ind = round(x);
            y_ind = round(y);

            image_cords(y_count,x_count,1) = y_ind;
            image_cords(y_count,x_count,2) = x_ind;

            x_count = x_count + 1;
        end
    end
    
    y_count = 0;

    % even offset by gap
    for y_o = h_start_even : h_advan*2 : im_in_h
        y_count = y_count + 2;
        x_count = 1;

        w_start_even = w_start_even_o + y_count*x_shift_const;
        for x = w_start_even : w_advan : im_in_w
            y = y_o + x_count*y_shift_const;
            x_ind = round(x);
            y_ind = round(y);

            image_cords(y_count,x_count,1) = y_ind;
            image_cords(y_count,x_count,2) = x_ind;

            x_count = x_count + 1;
        end
    end
    
    image_cords = image_cords(4:373, 4:323, :);
    img = Lytro_RAW; img(sub2ind(size(img), image_cords(:, :, 1), image_cords(:, :, 2))) = 1; figure;imshow(img)