import time
import argparse
import sys
import os
import json
import Network
import Dataset
import Train
import torch
import shutil

parser = argparse.ArgumentParser(
    description='''Python script for training a neural network to emulate time-varying effects that controlled
    by some kind of periodic LFO''')

# arguments for the training/test data locations and file names
parser.add_argument('--pedal', '-p', default='DonnerFlanger', help='Pedal to be modelled')

args = parser.parse_args()

if __name__ == "__main__":
    """The main method creates the LSTM network, trains it and carries out validation """

    with open('Results/FDonnerFlangerra12c12rg9Singles1106/config.json', 'r') as f:
        configs = json.load(f)
        for hypar in configs:
            args.__setattr__(hypar, configs[hypar])
    dirPath = 'Results/' + args.pedal + ''.join([s for s in args.load_config if s.isdigit()])

    # Create dataset object
    data = Dataset.DataSet(args)

    # Create instance of Network.RNN class
    network = Network.RNN(hidden_size=args.hidden_size, num_layers=args.num_layers, unit_type=args.unit_type,
                          input_size=args.in_size)

    network.load_state_dict(torch.load(dirPath + '/modelBest.pt', map_location='cpu'))

    Train.test(data, network, dirPath, args)