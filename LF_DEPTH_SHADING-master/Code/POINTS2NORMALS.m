function [normals] = POINTS2NORMALS(points)

normals = zeros(size(points));

for i = 1:size(points,1)
    for j = 1:size(points,2)
        point1 = points(i,j,:);
        if (i==size(points,1))
            point2 = points(i,j,:);
        else
            point2 = points(i+1,j,:);
        end
        if (j==size(points,2))
            point3 = points(i,j,:);
        else
            point3 = points(i,j+1,:);
        end
        normal = cross((point2 - point1),(point3 - point1));
        normal = normal / sqrt((normal(1)^2)+(normal(2)^2)+(normal(3)^2));
        normals(i,j,:) = normal;
    end
end

end