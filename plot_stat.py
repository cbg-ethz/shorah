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


def plot_freq_mp(files, threshold=0.5):
    ''' Plot frequency for files freq_file.csv produced by diri_sampler
    '''
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    import numpy as np
    import sys
    
    # match filename to set xticks
    import re
    rule = re.compile('w(\d*)-(\d*)')
    xt = []
    xv = []
    
    # Set figure scale, font and similia
    fig = plt.gca()
    fig.set_yscale('log', basey=2)
    fig.hold(True)
    
    for fnumber, filename in enumerate(files):
        print >> sys.stderr, 'doing ',  filename, 'number %d of %d' % (fnumber+1, len(files))
        with open(filename) as f:
            r = f.readline().split()[2:]
            history = len(r)
            dtype = [('haplotypes', 'S16')] + [('', np.int32)]*(history+1)
            y = np.loadtxt(f, dtype=dtype)
            y = y.view(np.dtype([('haplotypes', 'S16'), ('support', np.int32, 1),
                                 ('reads', np.int32, history)]))
        
        haplotypes = y['haplotypes']
        support = y['support']
        reads = y['reads']
    
        assert history == len(reads[0]), 'Some iteration is missing'
        
        cov = np.sum(reads, axis=0)
        coverage = cov[0]
        chk = np.equal(cov/cov[0], np.ones(history))
        if not np.all(chk):
            print 'Some read is missing in %s' % filename
            continue
        n_supported=0
        for i, h in enumerate(haplotypes):
            if support[i] >= threshold*history:
                n_supported = i
                
        reads_t = reads[:n_supported,:].transpose().astype(float)/coverage
        pos = [fnumber+1]*n_supported
        try:
            d = plt.boxplot(reads_t, positions=pos, sym='', hold=True)
        except:
            continue
        try:
            v = rule.search(filename)
            xt.append(v.group(1))
        except:
            xt.append('')
        xv.append(fnumber+1)
    plt.title('Boxplot for frequency of haplotypes')
    plt.xticks(xv, xt)
    fig = plt.gca()
    plt.xlim(0.5, fnumber+1.5)
    
    # save the figure
    imtype = 'pdf'
    plt.savefig('box_freq.%s' % imtype, dpi=None, facecolor='w', edgecolor='w',\
                    orientation='portrait', papertype=None, format=imtype,\
                    transparent=False)


def stat_freq(files, threshold=0.5):
    ''' If matplotlib is not available, use csv module to print some information
    '''
    import csv
    import sys
    from math import sqrt
    
    csv.register_dialect('support', delimiter='\t', skipinitialspace=True)
    on = open('freq_stat.dat', 'w')
    for filename in files:
        data = csv.reader(open(filename), dialect='support')
        support = []
        freq = []
        labels = []
        data.next()
        for row in data:
            labels.append(row[0])
            support.append(int(row[1]))
            try:
                freq.append(map(int, row[2:]))
            except:
                sys.exit(row)
        coverage = sum([f[0] for f in freq])
        history = len(freq[0])
        reads = {}
        std = {}
        
        on.write('#input file=%s, coverage=%d, haplotypes with > %d %% support:\n' % (filename, coverage, int(threshold*100)))
        for i, sup in enumerate(support):
            if float(support[i])/history > threshold:
                reads[labels[i]] = float(sum(freq[i]))/history
                std[labels[i]] = sqrt(float(sum([ (d-reads[labels[i]])**2 for d in freq[i] ]))/(history-1))
                on.write('%s has support %4.2f %% and frequency %5.2E  +/- %5.2E \n'%
                         (labels[i], 100*float(support[i])/history, reads[labels[i]]/coverage, std[labels[i]]/coverage))
        on.write('\n')
    on.close()


def plot_supp_mp(files):
    ''' Plot support only for files freq_file.csv produced by diri_sampler
    '''
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    import numpy as np
    import sys
    # match filename to set xticks
    import re
    rule = re.compile('w(\d*)-(\d*)')
    xt = []
    xv = []
    
    # Set figure scale, font and similia
    fig = plt.gca()
    fig.hold(True)
    plt.title('Top supported haplotypes')
    for fnumber, filename in enumerate(files):
        print >> sys.stderr, 'doing file ', filename, 'number %d of %d' % (fnumber+1, len(files))
        with open(filename) as f:
            r = f.readline().split()[2:]
            history = len(r)
            dtype = [('haplotypes', 'S16')] + [('', np.int32)]*(history+1)
            y = np.loadtxt(f, dtype=dtype)
            y = y.view(np.dtype([('haplotypes', 'S16'), ('support', np.int32, 1),
                                 ('reads', np.int32, history)]))

        haplotypes = y['haplotypes']
        support = y['support']
        normalized = support.astype(float)/history
        to_plot = []
        for t in normalized:
            if t > 0.1:
                to_plot.append(t)
        plt.plot(to_plot, 'o-', label=filename)
    plt.xlim(xmin=-0.2)
    plt.ylim(ymax=1)
    plt.xlabel('haplotype index')
    plt.ylabel('support')
    # plt.legend(loc=3, label='files')
    # save the figure
    imtype = 'pdf'
    plt.savefig('support.%s' % imtype, dpi=None, facecolor='w', edgecolor='w',\
                    orientation='portrait', papertype=None, format=imtype,\
                    transparent=False)


def main():
    ''' What does a main do?
    '''
    import sys
    import glob
    args = sys.argv
    
    try:
        directory = args[1]
        mode = args[2]
    except:
        sys.exit('\tusage: plot_stat.py directory mode[freq|support]')

    files = glob.glob('*-freq.csv')    
    if mode == 'freq':
        
        #try:
        import matplotlib
        plot_freq_mp(files)
        #except:
         #   stat_freq(files)
    
    if mode == 'support':
        plot_supp_mp(files)
        sys.exit()
        try:
            import matplotlib
            plot_supp_mp(files)
        except:
            sys.exit('try with mode freq')
    
if __name__ == "__main__":
    main()
