import os
import yaml


def remove_files():
    config = yaml.safe_load(open("config.yml"))
    if os.path.isfile(config['coder-generator']):
        os.remove(config['coder-generator'])
    if os.path.isfile(config['decoder-parity-check']):
        os.remove(config['decoder-parity-check'])
    if os.path.isfile(config['decoder-syndrome-decoding']):
        os.remove(config['decoder-syndrome-decoding'])


def read_file_to_list(name):
    lines = []
    with open(name) as file:
        for line in file:
            line = line.strip()
            lines.append(int(line, 2))
    return lines


def read_file_to_dict(name):
    dictionary = {}
    with open(name) as file:
        key = 0
        for line in file:
            line = line.rstrip()
            if line.isdigit():
                dictionary[key].append(int(line, 2))
            if line.startswith('S'):
                key = int(line.partition(':')[2], 2)
                dictionary[key] = []
    return dictionary


def print_files(g_matrix, h_matrix, syndrome_table, n, k):
    config = yaml.safe_load(open("config.yml"))
    r = n - k
    with open(config['Code_G_Matrix'], 'a') as file:
        [(file.write('{0:0>{width}b}\\n'.format(i, width=n)) for i in g_matrix)]
    with open(config['Decoder_H_Matrix'], 'a') as file:
        [(file.write('{0:0>{width}b}\\n'.format(i, width=n)) for i in h_matrix)]
    with open(config['Syndrome_Decoding_Table'], 'a') as file:
        for key in syndrome_table.keys():
            file.write('\n')
            file.write('S:{0:0>{width}b}'.format(key, width=r))
            file.write('\n')
            file.write('----------------Keys section-------------------')
            file.write('\n')
            for word in syndrome_table[key]:
                file.write('{0:0>{width}b}'.format(word, width=n))
                file.write('\n')
