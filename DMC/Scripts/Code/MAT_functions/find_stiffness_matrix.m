function K=find_stiffness_matrix(elements,coordinates)

K = zeros(length(coordinates));

% 1. Find local basis functions
for d=1:2  % x and y
    for i=1:3   % 3 local basis functions for every element
        % For each point in each element extract coordinates
        coord(d,i,:)=coordinates(elements(:,i),d);
    end
end   

K11 = 0.5*[1 -1 0; -1 1 0; 0 0 0];
K22 = 0.5*[1 0 -1; 0 0 0;-1 0 1];
K12 = 0.5*[1 0 -1; -1 0 1; 0 0 0];
for i = 1:length(coord)
    tri_coords = coord(:,:,i);
    B = [tri_coords(1,2) - tri_coords(1,1), tri_coords(1,3) - tri_coords(1,1); ...
         tri_coords(2,2) - tri_coords(2,1), tri_coords(2,3) - tri_coords(2,1)];
    C = inv(B)*(inv(B)');
    Ktri = abs(det(B))*(C(1,1)*K11 + C(2,2)*K22 + C(1,2)*(K12 + K12'));
    point_numbers = elements(i,:);
    for j = 1:3
        for k = 1:3
            K(point_numbers(j),point_numbers(k)) =  K(point_numbers(j),point_numbers(k)) + Ktri(j,k);
        end
    end
end

K = sparse(K);
