function [shading_normals] = SHADING2NORMALS(shading, lighting, LF_parameters)

%Parameters
y_size = LF_parameters.y_size;
x_size = LF_parameters.x_size;

M = [0.429043*lighting(9) 0.429043*lighting(5) 0.429043*lighting(8) 0.511664*lighting(4); 
     0.429043*lighting(5) -0.429043*lighting(9) 0.429043*lighting(6) 0.511664*lighting(2);
     0.429043*lighting(8) 0.429043*lighting(6) 0.743125*lighting(7) 0.511664*lighting(3);
     0.511664*lighting(4) 0.511664*lighting(2) 0.511664*lighting(3) ((0.886227*lighting(1)) - (0.247708*lighting(7)))];

shading_normals = zeros(y_size, x_size, 3);

options = optimoptions(@fminunc, 'Display', 'off');

for j = 1:y_size
    for i = 1:x_size
        %point_normal = lsqnonlin(@(x)((x'*M*x)-shading(j,i)),zeros(4,1),[],[],options);
        %point_normal = quadprog(M, [0 0 0 -shading(j,i)], [], [], [], [], [], [], [], options);
        point_normal = fminunc(@(x)((x'*M*x)-(256*shading(j,i)))^2,zeros(4,1),options);
        shading_normals(j,i,:) = point_normal(1:3);
    end
end

end