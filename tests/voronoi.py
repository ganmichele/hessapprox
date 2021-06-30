#!/usr/bin/env python
"""
Script to generate Voronoi map and NGas movie given input data
Michele Gandolfi 2021
"""

import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
import os
import sys
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
try:
    from hessapprox.Ngas import ngas
except Exception as e1:
    print( 'Assuming module is in ../hessapprox/. Trying to add to sys.path..')
    sys.path.append( '../')
    try:
        from hessapprox.Ngas import ngas
    except Exception as e:
        quit( 'Cannot import hessapprox package, Aborting..')


def update(frame):
    ln.set_data( frame[:,0], frame[:,1])
    if np.all( frame == ws[-1]):
        vor = Voronoi( frame)
        voronoi_plot_2d( vor
                        ,ax=ax
                        ,show_vertices=False
                        )
    return ln,

def init():
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ax.plot( t[::1,0], t[::1,1], '-ko', markersize=2.5)
    return ln,

if __name__ == '__main__':
    t = np.genfromtxt( 'nm_water/h2o_NMq.dat')
    t = t[:150,:2]
    t = (t - t.min(0)) / (t.max(0) - t.min(0))
    #n = np.random.rand( len( t)//10, 2)
    #n = t[::10,:]
    if os.path.isfile( 'weights29.dat'):
        w = np.genfromtxt( 'weights29.dat')
    else:
        ngas = ngas( nn=25, print_weights_epochs=True)
        ngas.train( t, init_method='ran', epochs=30)
        w = ngas.weights
    vor = Voronoi( w)

    fig = voronoi_plot_2d( vor, show_vertices=False)
    xlim = plt.xlim()
    ylim = plt.ylim()
    del fig

    fig, ax = plt.subplots()
    ax.set_xlim( xlim)
    ax.set_ylim( ylim)

    ln, = plt.plot([], [], 'ro',
                   markersize=15,
                   markeredgewidth=1.0,
                   markeredgecolor='k' )

    metadata = dict(title='NGas movie',
                    artist='Michele',
                    comment='NGas test movie')

    ws = np.zeros( (30, 25, 2))
    for i in range( 30):
        ws[i,:,:] = np.genfromtxt( 'weights{0}.dat'.format( i))

    print( 'Preparing animation..')
    ani = FuncAnimation( fig
                        ,update
                        ,frames=ws
                        ,init_func=init
                        ,fargs=None
                        ,blit=True)

    print( 'Writing movie..')
    ani.save('water2d_NG.mp4'
              ,fps=3
              ,dpi=500
              ,bitrate=-1
              ,extra_args=['-vcodec', 'libx264'])

    fig = voronoi_plot_2d( vor
              ,ax=ax
              ,point_size=25.0
              ,show_vertices=False)
    for i in range( 30):
        os.remove( 'weights{0}.dat'.format( i))
    print( 'Movie written in "water2d_NG.mp4"')
