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

    '''
    dataname = []
    Singles = 6
    for each in range(Singles):
        dataname.append('FDonnerFlangerra12c12rg12Singles' + str(each + 1))
    for each in range(Singles):
        dataname.append('FDonnerFlangerra12c12rg9Singles' + str(each + 1))
    '''

    '''
    dataname = []
    Singles = 6
    for each in range(Singles):
        dataname.append('FDonnerFlangerra12c12rg12Singles' + str(each + 1))
    for each in range(Singles):
        dataname.append('FDonnerFlangerra12c12rg9Singles' + str(each + 1))
    '''
    '''
    dataname = []
    spacing = [5,10,20,30,40,50]
    for each in spacing:
        dataname.append('SweetTonePhaser_spc' + str(each) + '_notch2')
    for each in spacing:
        dataname.append('SweetTonePhaser_spc' + str(each) + '_notch5')
    '''
    '''
    dataname = []
    Singles = 6
    for each in range(Singles):
        dataname.append('GDonnerFlangerra12c12rg12Singles' + str(each + 1))
    #for each in range(Singles):
    #    dataname.append('GDonnerFlangerra12c12rg9Singles' + str(each + 1))
    #for each in range(Singles):
    #    dataname.append('GDonnerFlangerra12c12rg3Singles' + str(each + 1))

    hids = [32, 48, 64]
    '''
    '''
    dataname = []
    Singles = 6
    for each in range(Singles):
        dataname.append('HDonnerFlangerra12c9rg9Singles' + str(each + 1))
    for each in range(Singles):
        for leach in range(each+1, Singles):
            dataname.append('HDonnerFlangerra12c9rg9Doubles' + str(each + 1) + str(leach + 1))
    #for each in range(Singles):
    #    dataname.append('GDonnerFlangerra12c12rg9Singles' + str(each + 1))
    #for each in range(Singles):
    #    dataname.append('GDonnerFlangerra12c12rg3Singles' + str(each + 1))
    
    '''
    '''
    dataname = []
    Singles = 3
    for each in range(Singles):
        dataname.append('KDonnerFlangerra12c12rg9Singles' + str(each + 1))
    for each in range(Singles):
        dataname.append('KDonnerFlangerra12c12rg12Singles' + str(each + 1))
    for each in range(Singles):
        dataname.append('KDonnerFlangerra12c12rg3Singles' + str(each + 1))

    hids = [32, 48]
    '''

    hids = [32, 48, 64]

    dataname = ['JDonnerFlangerDoubles_rg3_rg9','JDonnerFlangerDoubles_rg3_rg12','JDonnerFlangerDoubles_rg9_rg12',
                'JDonnerFlangerTriple', 'KDonnerFlangerDoubles_rg3_rg9','KDonnerFlangerDoubles_rg3_rg12',
                'KDonnerFlangerDoubles_rg9_rg12', 'KDonnerFlangerTriple']

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    i = 600

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