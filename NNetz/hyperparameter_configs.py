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

	hids = [4,8,16,32,64]	
	lays = [1::4]	

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    i = 0

    for hs in hids:
    	for lay in lays:

        	configs = {'batch_size': 8,
        				'hidden_size': hs,
        				'num_lays': lay,
                   		'epochs': 200}
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
