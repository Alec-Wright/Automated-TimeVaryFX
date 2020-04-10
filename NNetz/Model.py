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
parser.add_argument('--pedal', '-p', default='SwTo', help='Pedal to be modelled')
parser.add_argument('--data_location', '-dl', default='../Dataset', help='Location of the "Data" directory')
parser.add_argument('--config_location', '-cl', default='Configs', help='Location of the "Configs" directory')

# arguments for the processing of the training/test data
parser.add_argument('--seg_len', '-sl', type=int,   default=44100, help='Length of audio segments in samples')

# arguments for the network training
parser.add_argument('--learn_rate',   '-lr',  type=float, default=0.0005,  help='learning rate')
parser.add_argument('--epochs',       '-eps', type=int,   default=750,   help='Number of training epochs to run')
parser.add_argument('--batch_size',   '-bs',  type=int,   default=10,    help='Training mini-batch size')
parser.add_argument('--val_freq', '-vfr', type=int,   default=1,    help='Validation Frequency (in epochs)')
parser.add_argument('--loss_fcn',     '-lf',  default='ESR',        help='MSE or ESR or ESRPre')
parser.add_argument('--up_fr',  '-uf',  type=int,   default=2048, help='Number of samples to process between '
                                                                             'parameter updates')
parser.add_argument('--init_len',     '-il',  type=int,   default=1000, help='Number of sequence samples to process'
                                                                             'before starting weight updates ')
parser.add_argument('--val_chunk',     '-vs',  type=int,   default=200000, help='Number of sequence samples to process'
                                                                               'in each chunk of validation ')
parser.add_argument('--test_chunk',   '-tc',  type=int,   default=200000, help='Number of sequence samples to process'
                                                                               'in each chunk of validation ')
parser.add_argument('--pre_filt',   '-pc',   default=[1, -0.85], help='FIR filter coefficients to be used during training, '
                                                                  'set to [] for no pre-emph filter, takes csv file'
                                                                  'or set to "adap" to do adaptive pre-emph filtering')
parser.add_argument('--low_pass',   '-lp', type=int,  default=0, help='flag to determine is lowpass filtering is applied'
                                                              'in the loss function [1, 0.85] FIR filter')
parser.add_argument('--segments',   '-sgs', type=int,  default=[3])

# arguments for loading, if restarting training
parser.add_argument('--cur_epoch',   '-cue',   type=int,    default=0,  help='epoch number to start on')
parser.add_argument('--vloss_list', '-vlo',  type=list,   default=[], help='val losses so far')
parser.add_argument('--tloss_list', '-tlo',  type=list,   default=[], help='train losses so far')

# arguments for the network structure
parser.add_argument('--in_size',  '-is', default=2, type=int, help='Number of input dimensions')
parser.add_argument('--num_layers',  '-nl', default=1, type=int, help='Number of recurrent layers')
parser.add_argument('--hidden_size', '-hs', default=16, type=int, help='Recurrent hidden state size')
parser.add_argument('--unit_type',   '-ut', default='LSTM',      help='LSTM or GRU')

parser.add_argument('--load_config', '-l', type=str, help="Config file path, to a JSON config file, any arguments"
                                                          "listed in the config file will replace the default values"
                    , default=None)
parser.add_argument('--printing', '-pr', action="store", dest="prin", type=str,
                    help="set to true to turn on console output", default=True)

args = parser.parse_args()

if __name__ == "__main__":
    """The main method creates the LSTM network, trains it and carries out validation """

    startTime = time.time()

    if args.load_config:
        configPath = args.config_location + '/config' + str(args.load_config) + '.json'
        with open(configPath, 'r') as f:
            configs = json.load(f)
            for hypar in configs:
                args.__setattr__(hypar, configs[hypar])
        dirPath = 'Results/' + args.pedal + ''.join([s for s in args.load_config if s.isdigit()])

        if os.path.isfile(dirPath + '/config.json'):
            configPath = dirPath + '/config.json'
            with open(configPath, 'r') as f:
                configs = json.load(f)
                for hypar in configs:
                    args.__setattr__(hypar, configs[hypar])
    else:
        dirPath = 'Results/' + args.pedal + '_hs-' + str(args.hidden_size) + \
                  '_nl-' + str(args.num_layers)

    # Check if a cuda device is available
    if not torch.cuda.is_available():
        print('cuda device not available')
        cuda = 0
    else:
        torch.set_default_tensor_type('torch.cuda.FloatTensor')
        torch.cuda.set_device(0)
        print('cuda device available')
        cuda = 1

    # Create dataset object
    data = Dataset.DataSet(args)

    # Create instance of Network.RNN class
    network = Network.RNN(hidden_size=args.hidden_size, num_layers=args.num_layers, unit_type=args.unit_type,
                          input_size=args.in_size)

    # If training is restarting, load model parameters
    if args.cur_epoch:
        try:
            network.load_state_dict(torch.load(dirPath + '/model.pt', map_location='cpu'))
        except:
            print('existing network not found, starting training from epoch 0')
            args.cur_epoch = 0
    if args.cur_epoch == args.epochs:
        print('training already complete - exiting')
        #exit()

    # Otherwise create directory for output
    if not args.cur_epoch:
        try:
            os.mkdir(dirPath)
        except FileExistsError:
            shutil.rmtree(dirPath)
            os.mkdir(dirPath)

    if cuda:
        network = network.cuda()

    optimizer = torch.optim.Adam(network.parameters(), lr=args.learn_rate)

    Train.train(data, network, optimizer, dirPath, args)

    trainTime = time.time() - startTime
    with open(dirPath + "/train_time.txt", "w") as output:
        output.write(str(trainTime / 60) + ' mins')

    network.load_state_dict(torch.load(dirPath + '/modelBest.pt', map_location='cpu'))

    print('Running Test Set')

    Train.test(data, network, dirPath, args)