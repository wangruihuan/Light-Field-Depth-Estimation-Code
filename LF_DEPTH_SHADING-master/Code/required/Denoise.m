function im_out = Denoise(im, radius, color_sigma)

im_in_p   = (max(min(round(im*255),255),0))                    ;
im_in_p_d = bilateral_filter_c(im_in_p,radius,color_sigma)          ;
im_in_p   = im_in_p_d/255                                           ;
im_out   = imfilter(im_in_p,fspecial('unsharp'))                   ;

end