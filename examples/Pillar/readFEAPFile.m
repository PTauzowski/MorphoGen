function [nodes, elements] = readFEAPFile(filename)
    % Reads nodes and elements from a FEAP input file.
    %
    % Args:
    %    filename (string): Path to the FEAP file.
    %
    % Returns:
    %    nodes (matrix): Nodal coordinates matrix [node_id, x, y, z].
    %    elements (cell array): Element connectivity {elem_id, node_ids}.

    fromFEAPnum = [1 9 2 12 21 10  4 11 3 17 25 18 23 27 24 20 26 19  5  13 6 16 22 14 8 15 7];
    
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end

    % Read header information
    task_name = fgetl(fid); % Read task name
    metadata_line = fgetl(fid); % Read metadata line
    metadata = sscanf(metadata_line, '%d');
    
    num_nodes = metadata(1);
    num_elements = metadata(2);
    dim = metadata(4);
    dofs_per_node = metadata(5);
    nodes_per_element = metadata(6);

    % Initialize storage
    nodes = zeros(num_nodes, dim); % Columns: [node_id, x, y, z]
    elements = cell(num_elements, 2); % Columns: {element_id, node_ids}

    % Read nodal coordinates
    while ~feof(fid)
        line = fgetl(fid);
        if contains(line, 'COORdinates')
            break;
        end
    end

    for i = 1:num_nodes
        trash = fscanf(fid, ' %s', 2);
        nodes(i, :) = fscanf(fid, ' %f', 3)';
    end

    % Read element connectivity
    while ~feof(fid)
        line = fgetl(fid);
        if contains(line, 'ELEMents')
            break;
        end
    end

    elem_idx = 1; % Initialize element index at 1
    elements=zeros(num_elements,27);
    for k=1:num_elements
            trash = fscanf(fid, ' %s', 3);
            element = fscanf(fid, ' %d', 27)';
            elements(k,:) = element(fromFEAPnum)';
    end

    fclose(fid);
end

