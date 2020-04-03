import numpy as np
import soundfile as sf
import torch
import math

# Class which reads the data from either np or wav files
class DataSet:
    def __init__(self, args):
        try:
            # Try and load data from numpy files
            inp_train = np.load(args.data_location + '/npData/train/' + args.pedal + '-input.npy')
            tgt_train = np.load(args.data_location + '/npData/train/' + args.pedal + '-target.npy')
            inp_val = np.load(args.data_location + '/npData/val/' + args.pedal + '-input.npy')
            tgt_val = np.load(args.data_location + '/npData/val/' + args.pedal + '-target.npy')
            inp_test = np.load(args.data_location + '/npData/test/' + args.pedal + '-input.npy')
            tgt_test = np.load(args.data_location + '/npData/test/' + args.pedal + '-target.npy')
        except FileNotFoundError:
            try:
                print('numpy files not found - loading from wav')
                # Try and load data from wav files, then save as numpy files in npData directory
                inp_train, samplerate = sf.read(
                    args.data_location + '/Data/train/' + args.pedal + '-input.wav')
                np.save(args.data_location + '/npData/train/' + args.pedal + '-input', inp_train)

                tgt_train, samplerate = sf.read(
                    args.data_location + '/Data/train/' + args.pedal + '-target.wav')
                np.save(args.data_location + '/npData/train/' + args.pedal + '-target', tgt_train)

                inp_val, samplerate = sf.read(
                    args.data_location + '/Data/val/' + args.pedal + '-input.wav')
                np.save(args.data_location + '/npData/val/' + args.pedal + '-input', inp_val)

                tgt_val, samplerate = sf.read(
                    args.data_location + '/Data/val/' + args.pedal + '-target.wav')
                np.save(args.data_location + '/npData/val/' + args.pedal + '-target', tgt_val)

                inp_test, samplerate = sf.read(
                    args.data_location + '/Data/test/' + args.pedal + '-input.wav')
                np.save(args.data_location + '/npData/test/' + args.pedal + '-input', inp_test)

                tgt_test, samplerate = sf.read(
                    args.data_location + '/Data/test/' + args.pedal + '-target.wav')
                np.save(args.data_location + '/npData/test/' + args.pedal + '-target', tgt_test)

            except FileNotFoundError:
                print(args.data_location + '/npData/train/' + args.pedal + '-input.npy: ''Dataset not found - exiting')
                exit()

        if args.segments:
            segments = args.segments
            total_length = len(inp_train)
            seg_len = int(total_length/6)
            segments[0] -= 1

            inp_trainnew = inp_train[segments[0] * seg_len:(segments[0] + 1)*seg_len, :]
            tgt_trainnew = tgt_train[segments[0] * seg_len:(segments[0] + 1) * seg_len]
            inp_valnew   = inp_val[segments[0] * seg_len:(segments[0] + 1) * seg_len, :]
            tgt_valnew   = tgt_val[segments[0] * seg_len:(segments[0] + 1) * seg_len]
            inp_testnew  = inp_test[segments[0] * seg_len:(segments[0] + 1)*seg_len, :]
            tgt_testnew  = tgt_test[segments[0] * seg_len:(segments[0] + 1) * seg_len]

            if len(segments) > 1:
                for segs in segments[1:]:
                    segs -= 1
                    inp_trainnew = np.concatenate((inp_trainnew,inp_train[segs * seg_len:(segs + 1)*seg_len, :]), 0)
                    tgt_trainnew = np.concatenate((tgt_trainnew, tgt_train[segs * seg_len:(segs + 1) * seg_len]), 0)

                    inp_valnew = np.concatenate((inp_valnew, inp_val[segs * seg_len:(segs + 1) * seg_len, :]), 0)
                    tgt_valnew = np.concatenate((tgt_valnew, tgt_val[segs * seg_len:(segs + 1) * seg_len]), 0)

                    inp_testnew = np.concatenate((inp_testnew, inp_test[segs * seg_len:(segs + 1) * seg_len, :]), 0)
                    tgt_testnew = np.concatenate((tgt_testnew, tgt_test[segs * seg_len:(segs + 1) * seg_len]), 0)

            inp_train = inp_trainnew
            tgt_train = tgt_trainnew
            inp_val = inp_valnew
            tgt_val = tgt_valnew
            inp_test = inp_testnew
            tgt_test = tgt_testnew

        # I ended up hard coding the sample rate as its not saved into the numpy file
        self.fs = 44100
        if inp_train.ndim == 1:
            inp_train = np.expand_dims(inp_train, 1)
            inp_test = np.expand_dims(inp_test, 1)
            inp_val = np.expand_dims(inp_val, 1)
        args.in_size = inp_train.shape[1]

        self.train_set = framify(inp_train, tgt_train, args.seg_len)

        self.val_set = framify(inp_val, tgt_val, inp_val.shape[0])

        self.test_set = framify(inp_test, tgt_test, inp_test.shape[0])


# converts continuous audio into frames, and creates a torch tensor from them
def framify(inp, target, frame_len):

    # Calculate the number of segments the training data will be split into
    seg_num = math.floor(inp.shape[0] / frame_len)

    # Find the number of input dimensions
    in_size = inp.shape[1]

    # Initialise training and validation set tensor matrices
    dataset = (torch.empty((frame_len, seg_num, in_size)),
               torch.empty((frame_len, seg_num, 1)))

    # Load the audio for the training set
    for i in range(seg_num):
        dataset[0][:, i, :] = torch.from_numpy(inp[i * frame_len:(i + 1) * frame_len])
        dataset[1][:, i, 0] = torch.from_numpy(target[i * frame_len:(i + 1) * frame_len])

    return dataset