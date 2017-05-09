function M=find_mass_matrix_surface(elements,coordinates)

M = sparse(length(coordinates),length(coordinates));

% 1. Find local basis functions
for d=1:3  % x and y and z
    for i=1:3   % 3 local basis functions for every element
        % For each point in each element extract coordinates
        coord(d,i,:)=coordinates(elements(:,i),d);
    end
end   


for i = 1:length(coord)
    tri_coords = coord(:,:,i);
    d1 = sqrt((tri_coords(1,1) - tri_coords(1,2))^2 + (tri_coords(2,1) - tri_coords(2,2))^2  + (tri_coords(3,1) - tri_coords(3,2))^2);
    d2 = sqrt((tri_coords(1,1) - tri_coords(1,3))^2 + (tri_coords(2,1) - tri_coords(2,3))^2  + (tri_coords(3,1) - tri_coords(3,3))^2);
    d3 = sqrt((tri_coords(1,2) - tri_coords(1,3))^2 + (tri_coords(2,2) - tri_coords(2,3))^2  + (tri_coords(3,2) - tri_coords(3,3))^2);
    
    tri_coords = zeros(2,3);
    tri_coords(:,2) = [d1,0];
    xx = (d2^2 - d3^2 + d1^2)/(2*d1);
    tri_coords(:,3) = [xx,sqrt(d2^2 - xx^2)];
        
    B = [tri_coords(1,2) - tri_coords(1,1), tri_coords(1,3) - tri_coords(1,1); ...
         tri_coords(2,2) - tri_coords(2,1), tri_coords(2,3) - tri_coords(2,1)];
    Mtri = abs(det(B))/24*(eye(3) + ones(3));
    point_numbers = elements(i,:);
    for j = 1:3
        for k = 1:3
            M(point_numbers(j),point_numbers(k)) =  M(point_numbers(j),point_numbers(k)) + Mtri(j,k);
        end
    end
end
