#!/usr/bin/env python

# Copyright 2007, 2008, 2009
# Niko Beerenwinkel,
# Nicholas Eriksson,
# Moritz Gerstung,
# Lukas Geyrhofer,
# Osvaldo Zagordi,
# ETH Zurich

# This file is part of ShoRAH.
# ShoRAH is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ShoRAH is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with ShoRAH.  If not, see <http://www.gnu.org/licenses/>.

def make_plot(filename, quantity1, quantity2=None):
    ''' Plot quantity for ds-out file
    '''
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    import numpy as np

    with open(filename) as f:
        r = f.readline().split()
        dtype = [('', np.int32)]*3 + [('',np.float32)]*2
        y = np.loadtxt(f, dtype=dtype)
        y = y.view(np.dtype([('iter', np.int32), ('K', np.int32), ('untouched', np.int32),
                             ('theta', np.float32), ('gamma', np.float32)]))
        '''
        iter = y['iter']
        K = y['K']
        untouch = y['untouch']
        theta = y['theta']
        gamma = y['gamma']
        ''' 
    # Set figure scale, font and similia
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(y[quantity1], 'bv-')
    ax1.set_xlabel('iterations')
    # Make the y-axis label and tick labels match the line color.
    q1l = quantity1
    if quantity1 == 'gamma':
        q1l = r'$\Gamma$'
    if quantity1 == 'theta':
        q1l = r'$\theta$'
    ax1.set_ylabel(q1l, color='b')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')

    if quantity2:
        ax2 = ax1.twinx()
        ax2.plot(y[quantity2], 'rx-')
        q2l = quantity2
        if quantity2 == 'gamma':
            q2l = r'$\Gamma$'
        if quantity2 == 'theta':
            q2l = r'$\theta$'
        ax2.set_ylabel(q2l, color='r')
        for tl in ax2.get_yticklabels():
            tl.set_color('r')

    # save the figure
    plt.title('Analysis of the run, file %s' % filename)
    imtype = 'pdf'
    if quantity2:
        plt.savefig('run_%s_%s_%s.%s' % (filename, quantity1, quantity2, imtype), dpi=None, facecolor='w', edgecolor='w',\
                        orientation='portrait', papertype=None, format=imtype,\
                        transparent=False)
    else:
        plt.savefig('run_%s_%s.%s' % (filename, quantity1, imtype), dpi=None, facecolor='w', edgecolor='w',\
                        orientation='portrait', papertype=None, format=imtype,\
                        transparent=False)

def main():
    ''' What does a main do?
    '''
    import sys
    
    args = sys.argv
    poss_quant = ['theta', 'gamma', 'K', 'untouched']
    try:
        filename = args[1]
        quantity1 = args[2]
    except:
        sys.exit('usage: plot_sampling.py file quantity1 [optional: quantity2] [%s]' % '|'.join(poss_quant))
        
    try:
        quantity2 = args[3]
        if quantity2 not in poss_quant:
            sys.exit('usage: plot_sampling.py file quantity1 [optional: quantity2] [%s]' % '|'.join(poss_quant))
    except:
        quantity2 = None
        
    if quantity1 not in poss_quant:
        sys.exit('usage: plot_sampling.py file quantity1 [optional: quantity2] [%s]' % '|'.join(poss_quant))
        
    make_plot(filename, quantity1, quantity2)
if __name__ == "__main__":
    main()
