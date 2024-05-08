 '''Puts the padded arrays (matlab data in a cube form) into two columns per cell'''
    column_list = []
    for i in range(len(padded_arrays[:-9])):#this should be 94 when running fully
        array_data = concatenated_data[0, i]
        if array_data.shape[1] >= 2:  # Check if array has at least 2 columns
            # Extract the first two columns
            first_column = array_data[:, 0]
            second_column = array_data[:, 1]
            # Combine the first and second columns into a single array
            combined_columns = np.column_stack((first_column, second_column)) 
            column_list.append(combined_columns)
            
    #Takes the column list, find where the spectra overlap and cut them both in between.
    for i, order in enumerate(column_list):
        try:
            next_order = column_list[i+1]#Defining the second order
            wave_max = order[-1,0]
            wave_min = next_order[0,0]
            wave_split = (wave_max - wave_min)/2 #Finds the middle of the overlap
            column_list[i] = order[order[:,0] < wave_min + wave_split]#Cuts the end 
            column_list[i+1] = next_order[next_order[:,0] > wave_max - wave_split]#Cuts the beginning 
        except:
            None

    # Concatenate all arrays into one long array
    long_array = np.concatenate(column_list)
