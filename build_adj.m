function adj_matrix = build_adj(arg, num_cells, network_type)

%for each of the tree, line, and star networks, can specify the number of
%cells and network type and create the proper matrix for further calculations

if strcmp(network_type,'tree')
    
    %for tree case, explicitly set the 2,4,8,16 cell adj matrix
    if (num_cells == 2)
        
        adj_matrix = [-arg arg; 1-arg -(1-arg)];
        
    elseif (num_cells == 4)
        
        adj_matrix = [-2*arg arg arg 0;
            1-arg -1 0 arg;
            1-arg 0 -(1-arg) 0;
            0 1-arg 0 -(1-arg)];
        
    elseif (num_cells == 8)
        
        adj_matrix = [-3*arg arg arg 0 arg 0 0 0;
            1-arg -(1+arg) 0 arg 0 arg 0 0;
            1-arg 0 -1 0 0 0 arg 0;
            0 1-arg 0 -1 0 0 0 arg;
            1-arg 0 0 0 -(1-arg) 0 0 0;
            0 1-arg 0 0 0 -(1-arg) 0 0;
            0 0 1-arg 0 0 0 -(1-arg) 0;
            0 0 0 1-arg 0 0 0 -(1-arg)];
        
    elseif (num_cells == 16)
        
        adj_matrix = [-4*arg arg arg 0 arg 0 0 0 arg 0 0 0 0 0 0 0;
            1-arg -(1+2*arg) 0 arg 0 arg 0 0 0 arg 0 0 0 0 0 0;
            1-arg 0 -(1+arg) 0 0 0 arg 0 0 0 arg 0 0 0 0 0;
            0 1-arg 0 -(1+arg) 0 0 0 arg 0 0 0 arg 0 0 0 0;
            1-arg 0 0 0 -1 0 0 0 0 0 0 0 arg 0 0 0;
            0 1-arg 0 0 0 -1 0 0 0 0 0 0 0 arg 0 0;
            0 0 1-arg 0 0 0 -1 0 0 0 0 0 0 0 arg 0;
            0 0 0 1-arg 0 0 0 -1 0 0 0 0 0 0 0 arg;
            1-arg 0 0 0 0 0 0 0 -(1-arg) 0 0 0 0 0 0 0;
            0 1-arg 0 0 0 0 0 0 0 -(1-arg) 0 0 0 0 0 0;
            0 0 1-arg 0 0 0 0 0 0 0 -(1-arg) 0 0 0 0 0;
            0 0 0 1-arg 0 0 0 0 0 0 0 -(1-arg) 0 0 0 0;
            0 0 0 0 1-arg 0 0 0 0 0 0 0 -(1-arg) 0 0 0;
            0 0 0 0 0 1-arg 0 0 0 0 0 0 0 -(1-arg) 0 0;
            0 0 0 0 0 0 1-arg 0 0 0 0 0 0 0 -(1-arg) 0;
            0 0 0 0 0 0 0 1-arg 0 0 0 0 0 0 0 -(1-arg)];
        
    else
        
        %tree needs to be power of 2 (less than 32)
        adj_matrix = arg;
        disp('Error: expected power of 2 input up to 16');
        
    end
    
elseif strcmp(network_type, 'line') %for the line, can automate network
    
    adj_matrix = zeros(num_cells);
    
    adj_matrix = adj_matrix - eye(num_cells);
    adj_matrix(1,1) = -arg;
    adj_matrix(end,end) = -(1-arg);
    
    for i = 1:(num_cells-1)
        adj_matrix(i+1,i) = 1-arg;
        adj_matrix(i,i+1) = arg;
    end
    
elseif strcmp(network_type, 'star') %for the star, can automate network
    
    adj_matrix = zeros(num_cells);
    
    adj_matrix = adj_matrix - (1-arg)*eye(num_cells);
    adj_matrix(1,1) = -(num_cells-1)*arg;
    adj_matrix(2:end,1) = (1-arg);
    adj_matrix(1,2:end) = arg;
    
else
    
    adj_matrix = arg;
    disp('Invalid network structure')
    
end
end