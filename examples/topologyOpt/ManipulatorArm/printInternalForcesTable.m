function  printInternalForcesTable(Fel)
    % Row and column descriptions
    rowNames = {
        'Axial Force Start', 'Shear Force Y Start', 'Shear Force Z Start', ...
        'Torsion Moment X Start', 'Moment Y Start', 'Torsion Moment Z Start', ...
        'Axial Force End', 'Shear Force Y End', 'Shear Force Z End', ...
        'Torsion Moment X End', 'Moment Y End', 'Moment Z End'
    };
    
    columnNames = {'Bar 1', 'Bar 2', 'Bar 3', 'Bar 4', 'Bar 5', 'Bar 6', 'Bar 7'};
    
    % Create table
    internalForcesTable = array2table(Fel, 'VariableNames', columnNames, 'RowNames', rowNames);
    % Add title as a property
    internalForcesTable.Properties.Description = 'Internal Forces for 3D Frame Finite Elements';
    
    % Display the table
    disp(internalForcesTable);
end

