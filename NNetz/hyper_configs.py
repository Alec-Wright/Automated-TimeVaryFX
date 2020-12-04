import json, os, sys

'''
Example of how to write JSON config files for a hyperparameter search.
After you've written the files you can call them with something like
python train.py --config_file=<path-to-config>

The json files can be loaded to a python dict using
with open(path_to_config_file, 'r') as f:
    configs = json.load(f)
'''

def main(output_dir):

    '''dataname = ['EDonnerFlanger']

    frcut = [2,3,4,6]
    for each in frcut:
        dataname.append('EDonnerFlangerfrcut' + str(each))

    MxErcut = [2,5,10,15]
    for each in MxErcut:
        dataname.append('EDonnerFlangerMxErcut' + str(each))

    AvErcut = [200,400,800, 1600, 2000]
    for each in AvErcut:
        dataname.append('EDonnerFlangerAvErcut' + str(each))'''
    dataname = []
    Singles = 16
    for each in range(16):
        dataname.append('EDonnerFlangerSingles' + str(each + 1))

    hids = [16]



    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    i = 70

    for pedals in dataname:
        for hid in hids:
            configs = {'pedal': pedals,
                       'hidden_size': hid}
            filename = 'config{:01}.json'.format(i)
            with open(os.path.join(output_dir, filename), 'w') as f:
                json.dump(configs, f) # Write to dict to JSON
            i += 1

if __name__ == '__main__':
    usage = 'usage: python hyperparameter_configs.py <output-dir>'
    try:
        output_dir = sys.argv[1]
        if sys.argv[1] == '--help' or sys.argv[1] == '-h':
            print(usage)
            sys.exit(0)
    except IndexError:
        print(usage)
        sys.exit(1)
    main(output_dir)