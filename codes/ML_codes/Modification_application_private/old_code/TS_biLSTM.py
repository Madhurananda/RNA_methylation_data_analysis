
import os,sys
import torch
import math

sys.path.append( os.path.dirname(__file__) )

from TS_global import *

class BiRNN(torch.nn.Module):
    def __init__(self, input_size, hidden_size, num_layers, dropout_rate=g_dropout_rate): #, num_classes, dropout_rate=g_dropout_rate):
        super(BiRNN, self).__init__()
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        self.lstm = torch.nn.LSTM(input_size, hidden_size, num_layers, batch_first=True, bidirectional=True) ; #, dropout=dropout_rate)
        #self.fc = torch.nn.Linear(hidden_size*2, num_classes)  # 2 for bidirection
        #self.dropout = torch.nn.Dropout(p=dropout_rate)

    def forward(self, x, device):
        # Set initial states
        h0 = torch.zeros(self.num_layers*2, x.size(0), self.hidden_size).to(device) # 2 for bidirection
        c0 = torch.zeros(self.num_layers*2, x.size(0), self.hidden_size).to(device)

        # Forward propagate LSTM
        out, _ = self.lstm(x, (h0, c0))  # out: tensor of shape (batch_size, seq_length, hidden_size*2)
        # Decode the hidden state of the last time step

        # Decode the hidden state of the last time step
        #out = self.fc(out[:, :, :])
        return out

class SignalPred (torch.nn.Module):
    def __init__ (self, signal_dis_size=dis_size, size_hidden = 512, size_layers = 3, dropout_rate=g_dropout_rate, ispretrain=False, pretrainmean=False, pretrainbp=False):
        super().__init__();

        self.signal_dis_size = signal_dis_size
        self.ispretrain = ispretrain
        self.pretrainmean = pretrainmean
        self.pretrainbp = pretrainbp

        self.signalPred = BiRNN( signal_dis_size, size_hidden, size_layers); #, dropout_rate)
        #if not ispretrain:
        #    for param in self.signalPred.features.parameters:
        #        param.requires_grad = False;
        #para_sum = 0;
        #for l_name, l_para in self.signalPred.named_parameters():
        #    para_sum += l_para.nelement(); #l_para.data.size(0);
        #print("Total parameters: {}".format( para_sum ))

        self.mask_linear = torch.nn.Linear( size_hidden*2, signal_dis_size )
        self.mask_linear1 = torch.nn.Linear( size_hidden*2, size_hidden//8)
        self.mask_linear1_2 = torch.nn.Linear( size_hidden//8, 1)
        self.mask_logsoftmax = torch.nn.LogSoftmax(dim=-1)
        self.mask_softmax = torch.nn.Softmax(dim=-1)
        #self.mask_dropout = torch.nn.Dropout(p=dropout_rate)
        self.mask_dropout1 = torch.nn.Dropout(p=dropout_rate)

        self.mask_linear1_bp = torch.nn.Linear( size_hidden*2, size_hidden//8)
        self.mask_dropout1_bp = torch.nn.Dropout(p=dropout_rate)
        self.mask_linear1_2_bp = torch.nn.Linear( size_hidden//8, 4)

        #self.pred_linear = torch.nn.Linear( size_hidden, 2 )
        self.pred_linear = torch.nn.Linear( size_hidden*2, size_hidden//8 )
        self.pred_linear1 = torch.nn.Linear( size_hidden//8, 2 )
        #self.pred_softmax = torch.nn.LogSoftmax(dim=-1)
        self.pred_softmax = torch.nn.Softmax(dim=-1)
        self.pred_dropout = torch.nn.Dropout(p=dropout_rate)
        #self.pred_dropout1 = torch.nn.Dropout(p=dropout_rate)

    def forward(self, x, device, maskedPos = None):
        x = self.signalPred(x, device);

        if self.pretrainbp:
            if maskedPos==None:
                print("Error!!!! In pretrain, maskedPos must NOT be NONE")
            x = self.mask_dropout1_bp(self.mask_linear1_bp(x));
            x = self.mask_linear1_2_bp(x)
            return x[torch.arange(x.size(0)), maskedPos, :]

        if self.ispretrain:
            if maskedPos==None:
                print("Error!!!! In pretrain, maskedPos must NOT be NONE")
            if self.pretrainmean:
                x = self.mask_dropout1(self.mask_linear1(x));
                x = self.mask_linear1_2(x)
                return x[torch.arange(x.size(0)), maskedPos, :]

            #x = self.mask_dropout(self.mask_linear(x));
            x = self.mask_linear(x)
            return (self.mask_softmax( x)[torch.arange(x.size(0)), maskedPos, :], self.mask_logsoftmax( x)[torch.arange(x.size(0)), maskedPos, :])
        else:
            x = x[:, x.size(1)//2, :].squeeze();
            x = self.pred_dropout(self.pred_linear(x));
            x = self.pred_linear1(x)
            return self.pred_softmax ( x )

            x = self.pred_dropout(self.pred_linear(x));
            #print("In SignalPred.forward", x.shape)
            #x = self.pred_dropout1(self.pred_linear1(x));
            x = self.pred_linear1(x)
            #print("In SignalPred.forward", x.shape)
            #return x[:, x.size(1)//2, :]
            return self.pred_softmax ( x[:, x.size(1)//2, :] )


if __name__=='__main__':
    #model = BiRNN(input_size=100, hidden_size=512, num_layers=3, num_classes=1)
    model = SignalPred()
    print(model)

    para_sum = 0;
    for l_name, l_para in model.named_parameters():
        para_sum += l_para.nelement();
    print("Total parameters: {}".format( para_sum ))


