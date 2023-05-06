import os, sys
import numpy as np
import math
from read_trace import *
import pandas
from pygnuplot import gnuplot
from optparse import OptionParser    

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-t", "--time", dest="timespan", default=None,
                  help="specify the time frame span eg. 100:200");

    (options, args) = parser.parse_args()
    
    fname = args[0]
    dump = read_trace(fname)
    path = dump[b'path']
    print(path.shape)
    
    disp = [0, path.shape[0]]
    if (options.timespan):
        disp = [int(x) for x in options.timespan.split(':')]
    
    num_states = path.shape[1]
    num_frames = disp[1] - disp[0]

    gp = gnuplot.Gnuplot()
    gp.set('term qt font "arial,15"')
    gp.unset('ytics')

    gp.set('yrange [-0.2:{}]'.format(num_states-0.8))
    gp.set('xrange [{}-0.2:{}-0.8]'.format(disp[0], disp[1]))

    trellis = np.asarray(path.shape[0] * [[i for i in range(num_states)]])
    df = pandas.DataFrame(trellis)
    
    cmd = ['u 1:{} w p pt 6 ps 4 notitle'.format(i) for i in range(2,num_states+2)]

    #label the metric
    metric = dump[b'metric']
    set_cmd = [['label "{:.2f}" at first {},{} center font "arial,12" front'.format(metric[i][j], i, trellis[i][j]) for i in range(disp[0], disp[1])] for j in range(num_states)]
    for c in set_cmd:
        gp.set(*c)

    # plot the historical path
    
    trace_cmd = []    
    for t in range(disp[0]+1, disp[1]):
        for s in range(num_states):
            trace_cmd += ['arrow from first {},{} to first {},{} nohead front'.format(t, s, t-1, path[t][s])]

    gp.set(*trace_cmd)
    #print(path)
    
    gp.plot_data(df, *cmd)
    input('press enter to quit')
    





    
