import CoreAudioML.miscfuncs as miscfuncs
import CoreAudioML.training as training
import CoreAudioML.dataset as dataset
import CoreAudioML.networks as networks
import torch
import json
from torch.utils.tensorboard import SummaryWriter
import argparse
import time
import os
import shutil
from scipy.io.wavfile import write

parser = argparse.ArgumentParser(
    description='''Python script for training a neural network to emulate time-varying effects that controlled
    by some kind of periodic LFO''')

# arguments for the training/test data locations and file names
parser.add_argument('--device', '-p', default='BehPhaser', help='This label describes what device is being modelled')
parser.add_argument('--data_location', '-dl', default='../LFOMeasurement/MeasurementOutputs',
                    help='Location of the "Data" directory')
parser.add_argument('--file_name', '-fn', default='BehPhaser_set1',
                  help='The filename of the wav file to be loaded as the input/target data, the script looks for files'
                       'with the filename and the extensions -input.wav and -target.wav ')
parser.add_argument('--config_location', '-cl', default='Configs', help='Location of the "Configs" directory')

# arguments for the processing of the training/test data
parser.add_argument('--seg_len', '-sl', type=int,   default=44100, help='Length of audio segments in samples')

# arguments for the network training
parser.add_argument('--learn_rate',   '-lr',  type=float, default=0.0005,  help='learning rate')
# number of epochs and validation
parser.add_argument('--epochs', '-eps', type=int, default=500, help='Max number of training epochs to run')
parser.add_argument('--validation_f', '-vfr', type=int, default=5, help='Validation Frequency (in epochs)')
# TO DO
parser.add_argument('--validation_p', '-vp', type=int, default=500,
                  help='How many validations without improvement before stopping training, None for no early stopping')

parser.add_argument('--batch_size',   '-bs',  type=int,   default=6,    help='Training mini-batch size')
parser.add_argument('--loss_fcns', '-lf', default={'ESR': 1},
                  help='Which loss functions, ESR, ESRPre, DC. Argument is a dictionary with each key representing a'
                       'loss function name and the corresponding value being the multiplication factor applied to that'
                       'loss function, used to control the contribution of each loss function to the overall loss ')
parser.add_argument('--up_fr',  '-uf',  type=int,   default=1024, help='Number of samples to process between '
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

    start_time = time.time()

    if args.load_config:
        configPath = args.config_location + '/config' + str(args.load_config) + '.json'
        with open(configPath, 'r') as f:
            configs = json.load(f)
            for hypar in configs:
                args.__setattr__(hypar, configs[hypar])
        save_path = 'Results/' + args.pedal + ''.join([s for s in args.load_config if s.isdigit()])

        if os.path.isfile(save_path + '/config.json'):
            configPath = save_path + '/config.json'
            with open(configPath, 'r') as f:
                configs = json.load(f)
                for hypar in configs:
                    args.__setattr__(hypar, configs[hypar])
    else:
        save_path = 'Results/' + args.device + '_hs-' + str(args.hidden_size) + \
                  '_nl-' + str(args.num_layers) + '_' + args.file_name

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
    dataset = dataset.DataSet(data_dir=args.data_location)
    dataset.create_subset('train', frame_len=args.seg_len)
    dataset.create_subset('val')
    dataset.create_subset('test')

    #dataset = dataset.DataSet('data')
    #dataset.load_file('train/BehPhaserToneoffSingles1', 'train')
    #dataset.load_file('val/BehPhaserToneoffSingles1', 'val')
    #dataset.load_file('test/BehPhaserToneoffSingles1', 'test')

    dataset.load_file(args.file_name, ['train', 'val', 'test'], [0.75, 0.125, 0.125])

    # Create instance of Network.RNN class
    network = networks.SimpleRNN(hidden_size=args.hidden_size, num_layers=args.num_layers, unit_type=args.unit_type,
                          input_size=args.in_size)


    # Otherwise create directory for output
    if not args.cur_epoch:
        try:
            os.mkdir(save_path)
        except FileExistsError:
            shutil.rmtree(save_path)
            os.mkdir(save_path)

    if cuda:
        network = network.cuda()

    optimiser = torch.optim.Adam(network.parameters(), lr=args.learn_rate)
    loss_functions = training.LossWrapper(args.loss_fcns, args.pre_filt)
    train_track = training.TrainTrack()
    writer = SummaryWriter(os.path.join('runs', save_path[8:]))

    # If training is restarting, this will ensure the previously elapsed training time is added to the total
    init_time = time.time() - start_time + train_track['total_time'] * 3600
    # Set network save_state flag to true, so when the save_model method is called the network weights are saved
    network.save_state = True
    patience_counter = 0

    # This is where training happens
    # the network records the last epoch number, so if training is restarted it will start at the correct epoch number
    for epoch in range(train_track['current_epoch'] + 1, args.epochs + 1):
        ep_st_time = time.time()

        # Run 1 epoch of training,
        epoch_loss = network.train_epoch(dataset.subsets['train'].data['input'][0],
                                         dataset.subsets['train'].data['target'][0],
                                         loss_functions, optimiser, args.batch_size)
        # Run validation
        if epoch % args.validation_f == 0:
            val_ep_st_time = time.time()
            val_output, val_loss = network.process_data(dataset.subsets['val'].data['input'][0],
                                                        dataset.subsets['val'].data['target'][0],
                                                        loss_functions, args.val_chunk)
            if val_loss < train_track['best_val_loss']:
                patience_counter = 0
                network.save_model('model_best', save_path)
                write(os.path.join(save_path, "best_val_out.wav"),
                      dataset.subsets['test'].fs, val_output.cpu().numpy()[:, 0, 0])
            else:
                patience_counter += 1
            train_track.val_epoch_update(val_loss.item(), val_ep_st_time, time.time())
            writer.add_scalar('Loss/val', train_track['validation_losses'][-1], epoch)

        train_track.train_epoch_update(epoch_loss.item(), ep_st_time, time.time(), init_time, epoch)
        # write loss to the tensorboard (just for recording purposes)
        writer.add_scalar('Loss/train', train_track['training_losses'][-1], epoch)

        network.save_model('model', save_path)
        miscfuncs.json_save(train_track, 'training_stats', save_path)

        if args.validation_p and patience_counter > args.validation_p:
            print('validation patience limit reached at epoch ' + str(epoch))
            break

    lossESR = training.ESRLoss()
    test_output, test_loss = network.process_data(dataset.subsets['test'].data['input'][0],
                                                  dataset.subsets['test'].data['target'][0], loss_functions,
                                                  args.test_chunk)
    test_loss_ESR = lossESR(test_output, dataset.subsets['test'].data['target'][0])
    write(os.path.join(save_path, "test_out_final.wav"), dataset.subsets['test'].fs, test_output.cpu().numpy()[:, 0, 0])
    writer.add_scalar('Loss/test_loss', test_loss.item(), 1)
    writer.add_scalar('Loss/test_lossESR', test_loss_ESR.item(), 1)
    train_track['test_loss_final'] = test_loss.item()
    train_track['test_lossESR_final'] = test_loss_ESR.item()

    best_val_net = miscfuncs.json_load('model_best', save_path)
    network = networks.load_model(best_val_net)
    test_output, test_loss = network.process_data(dataset.subsets['test'].data['input'][0],
                                                  dataset.subsets['test'].data['target'][0], loss_functions,
                                                  args.test_chunk)
    test_loss_ESR = lossESR(test_output, dataset.subsets['test'].data['target'][0])
    write(os.path.join(save_path, "test_out_bestv.wav"),
          dataset.subsets['test'].fs, test_output.cpu().numpy()[:, 0, 0])
    writer.add_scalar('Loss/test_loss', test_loss.item(), 2)
    writer.add_scalar('Loss/test_lossESR', test_loss_ESR.item(), 2)
    train_track['test_loss_best'] = test_loss.item()
    train_track['test_lossESR_best'] = test_loss_ESR.item()
    miscfuncs.json_save(train_track, 'training_stats', save_path)
    if cuda:
        with open(os.path.join(save_path, 'maxmemusage.txt'), 'w') as f:
            f.write(str(torch.cuda.max_memory_allocated()))

