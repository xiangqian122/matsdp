def read_proxigram_csv(proxigram_csv_dir):
    '''
    read the proxigram csv file.
    '''
    import numpy as np
    from .. import funcs
    with open(proxigram_csv_dir,'r') as f:
        line = f.readlines()
        line_num = len(line)
        row_num = len(funcs.split_line(line = line[0], separator = ','))
        data_set = np.array([None]*(line_num-1)*row_num)
        data_set.shape = (line_num-1),row_num
        for i_line in range(0,line_num-1):
            for i_row in range(0,row_num-1):
                data_set[i_line,i_row] = float(funcs.split_line(line = line[i_line+1], separator = ',')[i_row])
    elmt_num = int((row_num - 2) / 2)
    elmtname_list = []
    for i_elmt in range(0,elmt_num):
        elmtname_list.append(funcs.split_line(line = line[0], separator = ',')[i_elmt*2+1].split(' ')[0])       
    return data_set, elmtname_list
