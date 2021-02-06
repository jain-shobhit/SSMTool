function R = rotation_matrix(nodes)
% finds appropriate rotation matrix to rotate the elements which have
% dimension lower than the ambient dimension
% nodes: coordinates (in 1D, 2D or 3D) of the nodes in an element as
%           [x1 y1 z1;   node-1
%            x2 y2 z2;   node-2
%            ...   
%            xn yn zn];  node-n

dim = size(nodes,2);
elementdim = rank(diff(nodes));

switch dim
    case 1
        R = 1;
        
    case 2
        if elementdim == 1
            node1=nodes(1,:);
            node2=nodes(2,:);
            l=norm(node2-node1);
            e_x = (node2-node1)/l;
            R = [e_x(1) e_x(2);
                 -e_x(2) e_x(1)];
        else
            R = eye(2);
        end      
        
    case 3
        switch elementdim
            case 1
                node1 = nodes(1,:);
                node2 = nodes(2,:);
                d2  = node2-node1;
                l=norm(d2);                
                e_x = d2/l;                
                % find a third node (not aligned with node1 and node2
                % equation of plane normal to d2: d2(1)x + d2(2)y + d2(3)z = 0
                if d2(1)
                   n = [-(d2(2) + d2(3))/d2(1), 1 ,1];
                elseif d2(2)
                   n = [1, -(d2(1) + d2(3))/d2(2) ,1];                    
                else 
                    n = [1, 1, -(d2(1) + d2(2))/d2(3)];
                end                
                e_n = n/norm(n);
                node3 = node2 + l*e_n;
                
                d3  = node3-node1;
                e_z = cross(d2,d3)/norm(cross(d2,d3));
                e_y = cross(e_z,e_x);
                
                R = [ e_x(1) e_x(2) e_x(3);
                    e_y(1) e_y(2) e_y(3);
                    e_z(1) e_z(2) e_z(3)];               
                
            case 2
                node1=nodes(1,:);
                node2=nodes(2,:);
                node3=nodes(3,:);
                
                l=norm(node2-node1);
                
                e_x = (node2-node1)/l;
                d2  = node2-node1;
                d3  = node3-node1;
                e_z = cross(d2,d3)/norm(cross(d2,d3));
                e_y = cross(e_z,e_x);
                
                R = [ e_x(1) e_x(2) e_x(3);
                    e_y(1) e_y(2) e_y(3);
                    e_z(1) e_z(2) e_z(3)];               
            case 3
                R = eye(3);
        end
end

end