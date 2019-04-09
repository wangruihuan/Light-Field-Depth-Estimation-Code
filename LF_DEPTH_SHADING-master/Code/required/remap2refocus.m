directoryPath = 'analysis/spec_1/remap/';
Listing = dir(fullfile(directoryPath, '*.png'));
names = {Listing.name};
names = sort(names);
mkdir('analysis/spec_1/refocus');

for i = 1:256
    im = im2double(imread(strcat(directoryPath, names{i})));
    im_out = zeros(381, 330, 3);
    j = 1:381;
    k = 1:330;
    for m = 1:7;
        for n = 1:7
            im_out(j, k, :) = im_out(j, k, :) + im((j-1)*7+m, (k-1)*7+n, :);
        end
    end
    im_out = im_out/49;
    imwrite(im_out, strcat('analysis/spec_1/refocus/', names{i}));
end