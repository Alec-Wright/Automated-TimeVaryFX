import torch
import math
import numpy as np
import soundfile as sf
import json
import torch.nn as nn
import csv

def train(data, network, optimizer, dirPath, args):
    print('Starting training')

    # If no pre-filt is provided, the normal unfiltered ESR loss function is used
    if not args.pre_filt:
        loss_fnESR = ESRLoss()
    # If the pre-filt argument provided is a string, a csv file will be loaded to obtain the filter coefficients
    elif type(args.pre_filt) == str:
        with open('../' + args.pre_filt) as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            args.pre_filt = list(reader)
            args.pre_filt = args.pre_filt[0]
            for item in range(len(args.pre_filt)):
                args.pre_filt[item] = float(args.pre_filt[item])
        loss_fnESR = ESRLossPreEmph(args.pre_filt, args.low_pass)
    # If the pre-filt argument is simply filter coefficients, the network will use those coefficients
    else:
        loss_fnESR = ESRLossPreEmph(args.pre_filt, args.low_pass)

    loss_fnDC = DCLoss()

    if not args.vloss_list:
        best_val_loss = 100000
    else:
        best_val_loss = torch.tensor(min(args.vloss_list))

    for i in range(args.cur_epoch, args.epochs):
        print('Starting epoch ' + str(i+1))
        train_losses = []
        # shuffle the segments at the start of the epoch
        shuffle = torch.randperm(data.train_set[0].shape[1])

        for n in range(math.ceil(data.train_set[0].shape[1]/args.batch_size)):
            # Load batch of randomly selected segments
            batch = data.train_set[0][:, shuffle[n*args.batch_size:(n+1)*args.batch_size], :]
            target = data.train_set[1][:, shuffle[n*args.batch_size:(n+1)*args.batch_size], :]

            # Initialise the network hidden state then run initialisation samples through the network
            network.init_hidden(batch.shape[1])
            network(batch[0:args.init_len, :, :])

            # Zero the gradient buffers
            network.zero_grad()

            # Iterate over the remaining samples in the sequence batch
            for k in range(math.ceil((batch.shape[0]-args.init_len)/args.up_fr)):
                output = network(batch[args.init_len+(k*args.up_fr):args.init_len+((k+1)*args.up_fr), :, :])
                # Calculate loss
                loss = loss_fnESR(output, target[args.init_len+(k*args.up_fr):args.init_len+((k+1)*args.up_fr), :, :]) +\
                    loss_fnDC(output, target[args.init_len+(k*args.up_fr):args.init_len+((k+1)*args.up_fr), :, :])

                train_losses.append(loss.item())

                loss.backward()
                optimizer.step()

                # Set the network hidden state, to detach it from the computation graph
                network.set_hidden(network.hidden)
                network.zero_grad()
            #print(n/math.ceil(data.train_set[0].shape[1]/args.batch_size), 'loss:',loss)

        args.tloss_list.append(np.mean(train_losses))
        args.cur_epoch += 1

        # Run validation set
        if args.cur_epoch % args.val_freq == 0:
            print('Starting Validation')
            network.zero_grad()
            with torch.no_grad():
                # Initialise hidden state
                network.init_hidden(data.val_set[0].shape[1])
                #val_output = torch.empty(data.val_set[1].shape)
                val_output = network(data.val_set[0])
                # Process validation set
                #for l in range(int(val_output.size()[0]/args.val_chunk)):
                #    val_output[l*args.val_chunk:(l+1)*args.val_chunk] = network(data.val_set[0][l*args.val_chunk:(l+1)*args.val_chunk])
                #    network.set_hidden(network.hidden)
                # If the validation set doesn't divide evenly into the validation chunk length, process the remainder
                #if not (val_output.size()[0]/args.val_chunk).is_integer():
                #    val_output[(l+1)*args.val_chunk:-1] = network(data.val_set[0][(l+1)*args.val_chunk:-1])

                # Calculate the losses
                val_loss = loss_fnESR(val_output, data.val_set[1])
                val_loss += loss_fnDC(val_output, data.val_set[1])

                # Save the loss to the loss list
                args.vloss_list.append(val_loss.item())

            if val_loss < best_val_loss:
                best_val_loss = val_loss
                torch.save(network.state_dict(), dirPath + '/modelBest.pt')
                sf.write(
                    dirPath + '/bestvaloutput.wav', val_output.cpu().numpy()[:, 0, 0], 44100)

            print('epoch', i + 1, ', complete, validation loss:', val_loss.item(), 'best validation loss: ',
                  best_val_loss.item())

            torch.save(network.state_dict(), dirPath + '/model.pt')
            sf.write(
                dirPath + '/valoutput.wav', val_output.cpu().numpy()[:, 0, 0], 44100)

            np.savetxt(dirPath + "/vlosses.csv", args.vloss_list, delimiter=",")
            with open(dirPath + '/bestvloss.txt', 'w') as f:
                f.write(str(min(args.vloss_list)))
            #np.savetxt(dirPath + "/bestvloss.txt", max(args.vloss_list))
            np.savetxt(dirPath + "/trainlosses.csv", args.tloss_list, delimiter=",")

            with open(dirPath + "/config.json", "w") as output:
                output.write(json.dumps(vars(args)))


def test(data, network, dirPath, args):

    tst_ch = args.test_chunk
    network.init_hidden(args.cond_vals)
    filts = 0
    filtstest = 0

    with torch.no_grad():

        loss_fn = ESRLoss()

        if not args.pre_filt:
            loss_fnESR = ESRLoss()
        elif args.pre_filt == 'adap':
            print('adaptive filtering active')
            loss_fnESR = ESRLossAdapPreEmph()
            filtstest = data.test_set[2]
        elif type(args.pre_filt) == str:
            with open('Configs/' + args.pre_filt) as csvfile:
                reader = csv.reader(csvfile, delimiter=',')
                args.pre_filt = list(reader)
                args.pre_filt = args.pre_filt[0]
                for item in range(len(args.pre_filt)):
                    args.pre_filt[item] = float(args.pre_filt[item])
            loss_fnESR = ESRLossPreEmph(args.pre_filt, args.low_pass)
        else:
            loss_fnESR = ESRLossPreEmph(args.pre_filt, args.low_pass)

        loss_fnDC = DCLoss()

        output = torch.empty(data.test_set[1].shape)

        #network.load_state_dict(torch.load(dirPath + '/modelBest.pt'))

        for l in range(int(output.size()[0] / tst_ch)):
            output[l * tst_ch:(l + 1) * tst_ch] = network(data.test_set[0][l * tst_ch:(l + 1) * tst_ch])
            network.set_hidden(network.hidden)

        if not (output.size()[0] / tst_ch).is_integer():
            output[(l + 1) * tst_ch:-1] = network(data.test_set[0][(l + 1) * tst_ch:-1])

        test_loss = loss_fn(output, data.test_set[1])
        #test_loss_clf = loss_fnESR(output, data.test_set[1]) + loss_fnDC(output, data.test_set[1])

        test_lossclf = loss_fnESR(output, data.test_set[1])

        test_lossclf += loss_fnDC(output, data.test_set[1])

        if args.cond_vals > 1:
            for i in range(args.cond_vals):

                sf.write(dirPath + '/testoutput' + str(i) + '.wav', output.cpu().numpy()[:, i, 0], 44100)

            flat = output.cpu().numpy().flatten('F')
            sf.write(dirPath + '/outputcat.wav', flat, 44100)
        else:
            testlosses = np.empty([1])
            testlosses[0] = test_loss.item()
            testlossesclf = np.empty([1])
            testlossesclf[0] = test_lossclf.item()
            sf.write(dirPath + '/testoutput.wav', output.cpu().numpy()[:, 0, 0], 44100)
            with open(dirPath + '/tloss.txt', 'w') as f:
                f.write(str(testlosses))
            with open(dirPath + '/tlossclf.txt', 'w') as f:
                f.write(str(testlossesclf))

    with open(dirPath + "/config.json", "w") as output:
        output.write(json.dumps(vars(args)))

# Loss Functions
class ESRLoss(nn.Module):
    def __init__(self):
        super(ESRLoss, self).__init__()
        self.epsilon = 0.00001

    def forward(self, output, target):
        loss = torch.add(target, -output)
        loss = torch.pow(loss, 2)
        loss = torch.mean(loss)
        energy = torch.mean(torch.pow(target, 2)) + self.epsilon
        loss = torch.div(loss, energy)
        return loss


class DCLoss(nn.Module):
    def __init__(self):
        super(DCLoss, self).__init__()
        self.epsilon = 0.00001

    def forward(self, output, target):
        loss = torch.pow(torch.add(torch.mean(target, 0), -torch.mean(output, 0)), 2)
        loss = torch.mean(loss)
        energy = torch.mean(torch.pow(target, 2)) + self.epsilon
        loss = torch.div(loss, energy)
        return loss


# Pre emph loss includes a pre-emphasis filter
class ESRLossPreEmph(nn.Module):
    def __init__(self, preFilt, lp):
        super(ESRLossPreEmph, self).__init__()
        self.epsilon = 0.00001
        self.preFilt = preFilt
        self.zPad = len(preFilt) - 1
        self.lp = lp

    def forward(self, output, target):

        conFil = nn.Conv1d(1, 1, 2, bias=False)
        conFil.weight.data = torch.tensor([[self.preFilt]], requires_grad=False)

        output = torch.cat((torch.zeros(self.zPad, output.shape[1], 1), output))
        target = torch.cat((torch.zeros(self.zPad, target.shape[1], 1), target))

        output = conFil(output.permute(1, 2, 0))
        target = conFil(target.permute(1, 2, 0))

        if self.lp:
            lpFil = nn.Conv1d(1, 1, 2, bias=False)
            lpFil.weight.data = torch.tensor([[[0.85, 1]]], requires_grad=False)
            output = lpFil(output)
            target = lpFil(target)

        output = output.permute(2, 0, 1)
        target = target.permute(2, 0, 1)

        loss = torch.add(target, -output)
        loss = torch.pow(loss, 2)
        loss = torch.mean(loss)
        energy = torch.mean(torch.pow(target, 2)) + self.epsilon
        loss = torch.div(loss, energy)
        return loss