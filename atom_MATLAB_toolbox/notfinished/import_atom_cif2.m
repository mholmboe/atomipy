function cif_data = import_atom_cif2(filename)
    % Check if the file is a .cif file
    if ~endsWith(filename, '.cif')
        error('The file must be a .cif file');
    end

    % Open the CIF file
    fileID = fopen(filename, 'r');
    if fileID == -1
        error('File could not be opened.');
    end

    % Initialize a structure to store CIF data
    cif_data = struct;

    try
        % Read and parse the file line by line
        while ~feof(fileID)
            line = fgetl(fileID);

            % Example: Parse atom_site_fract_x, atom_site_fract_y, atom_site_fract_z
            if contains(line, 'atom_site_fract_x')
                cif_data.atom_site_fract_x = str2double(regexp(line, '\d+\.?\d*', 'match'));
            elseif contains(line, 'atom_site_fract_y')
                cif_data.atom_site_fract_y = str2double(regexp(line, '\d+\.?\d*', 'match'));
            elseif contains(line, 'atom_site_fract_z')
                cif_data.atom_site_fract_z = str2double(regexp(line, '\d+\.?\d*', 'match'));
            end

            % Add more parsing logic for other CIF properties as needed

        end
    catch
        fclose(fileID);
        error('Error occurred while reading the file.');
    end

    % Close the file
    fclose(fileID);
end