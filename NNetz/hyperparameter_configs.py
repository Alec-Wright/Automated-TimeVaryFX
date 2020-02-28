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

    pedal = ['HT5M10', 'Mesa550Crunch3']
    epochs = [50,100,200,500]

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    i = 0

    for peds in pedal:
    	for eps in epochs:

        	configs = {'pre_filt': 'b_Awght_mk2.csv',
            	        'batch_size': 8,
                	    'pedal': peds,
                   		'low_pass': 1,
                   		'epochs': eps}
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
