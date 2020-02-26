import torch
import torch.nn as nn

# Recurrent Neural Network class
class RNN(nn.Module):
    def __init__(self, input_size=1, hidden_size=16, output_size=1, num_layers=1, unit_type='LSTM'):
        super(RNN, self).__init__()
        # Set number of layers  and hidden_size for network layer/s
        self.input_size = input_size
        self.hidden_size = hidden_size
        self.output_size = output_size
        self.num_layers = num_layers
        self.unit_type = unit_type
        self.hidden = None

        self.layers = nn.ModuleList()
        # Create recurrent layer of the type specified by unit_type
        cell_types = {'LSTM': nn.LSTM, 'GRU': nn.GRU}
        try:
            for n in range(self.num_layers):
                self.layers.append(BasicBlock(self.hidden_size, self.input_size))
                input_size = output_size
        except KeyError:
            print('Invalid unit type, choose "LSTM" or "GRU"')
            exit()

    # Define forward pass, where x and y are the network inputs/outputs
    def forward(self, x):
        for n in range(self.num_layers):
            x = self.layers[n](x)
        return x

    # Initialise hidden state with zeros
    def init_hidden(self, batch_size):
        for n in range(self.num_layers):
            self.layers[n].init_hidden(batch_size)

    # Set hidden state to specified values, this resets gradient tracking on the hidden state
    def set_hidden(self, values):
        for n in range(self.num_layers):
            self.layers[n].set_hidden()


class BasicBlock(nn.Module):

    def __init__(self, hid_sz, input_size):
        super(BasicBlock, self).__init__()
        self.rec = nn.LSTM(input_size=input_size, hidden_size=hid_sz, num_layers=1)
        self.lin = nn.Linear(hid_sz, 1)
        self.unit_type = 'LSTM'
        self.hidden_size = hid_sz
        self.hidden = None

    def forward(self, x):
        res = x
        x, self.hidden = self.rec(x, self.hidden)
        x = self.lin(x) + res[:, :, 0:1]

        return x

    # Initialise hidden state with zeros
    def init_hidden(self, batch_size):
        self.hidden = torch.zeros(1, batch_size, self.hidden_size)
        if self.unit_type == 'LSTM':
            self.hidden = [self.hidden, torch.zeros(1, batch_size, self.hidden_size)]

    # Set hidden state to specified values, this resets gradient tracking on the hidden state
    def set_hidden(self):
        values = self.hidden
        if self.unit_type == 'LSTM':
            self.hidden = [values[0].clone().detach(), values[1].clone().detach()]
        elif self.unit_type == 'GRU':
            self.hidden = values.clone().detach()