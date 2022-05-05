def get_interface_positions(dataframe, column1, column2):
    interface_positions = {}
    for i in dataframe.index:
        if dataframe.at[i, column1] not in interface_positions and dataframe.at[i, column1 + '_IRES'] != '[]':
            interface_positions[dataframe.at[i, column1]] = dataframe.at[i, str(column1 + '_IRES')]
        elif dataframe.at[i, column1] in interface_positions and dataframe.at[i, column1 + '_IRES'] != '[]':
            interface_positions[dataframe.at[i, column1]] = interface_positions[dataframe.at[i, column1]].strip(
                ']') + ',' + (dataframe.at[i, str(column1 + '_IRES')]).strip('[')
        if dataframe.at[i, column2] not in interface_positions and dataframe.at[i, column2 + '_IRES'] != '[]':
            interface_positions[dataframe.at[i, column2]] = dataframe.at[i, str(column2 + '_IRES')]
        elif dataframe.at[i, column2] in interface_positions and dataframe.at[i, column2 + '_IRES'] != '[]':
            interface_positions[dataframe.at[i, column2]] = interface_positions[dataframe.at[i, column2]].strip(
                ']') + ',' + (dataframe.at[i, str(column2 + '_IRES')]).strip('[')

    try:
        for key, value in interface_positions.items():
            n = []
            m = []
            if value != '[]':
                valueList = value.split(',')
                valueList[0] = str(valueList[0]).strip('[')
                valueList[-1] = str(valueList[-1]).strip(']')
                for val in valueList:
                    if '-' in val:
                        for r in range(int(val.split('-')[0]), int(val.split('-')[1]) + 1):
                            n.append(r)
                    else:
                        m.append(int(val))
                fin = m + n

                interface_positions[key] = fin
    except:
        ValueError

    return interface_positions
