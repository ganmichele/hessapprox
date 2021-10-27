#!/usr/bin/env python
"""
Hessian Approximation Methods. Full notice in LICENSE file
Copyright (C) 2021 Michele Gandolfi, Michele Ceotto

NGas method as explained in
Gandolfi M, Ceotto M, J. Chem. Theory Comput. (2021)

Michele Gandolfi 2021
"""
import numpy as np
import sys
try:
    from calc_dist import calc_distance as fast_distance
except Exception as e:
    #print( 'Could not import routine for fast distance calculation')
    try:
        from scipy.spatial import distance_matrix as scipy_distance
    except Exception as e:
        print( 'Could not import any routine for fast distance calculation')


def pairwise_dist( a, b, dist_fun='eu'):
    """
    Compute pairwise distance between instances "a" and "b"
    """
    if dist_fun == 'cb':
        return abs( a - b).sum()
    else: # euclidean
        return np.sqrt( (a-b) @ (a-b) )

def batch_distance( a, b, dist_fun='eu'):
    """
    Compute distance between data in 'a' w.r.t data in 'b'
    """
    if 'calc_dist' in sys.modules: # Use faster Fortran code if imported
        if a.shape != b.shape:
            same = False
        else:
            same = np.array_equal( a, b)#(a==b).all()
        a.dtype, b.dtype = np.float64, np.float64
        a = np.asfortranarray( a) 
        b = np.asfortranarray( b) 
        return fast_distance( a, b, same, dist_fun)
    elif 'scipy' in sys.modules:
        p = 1 if dist_fun=='cb' else 2
        return scipy_distance( a, b, p=p)
    else:
        dist = np.zeros( ( len(a), len(b)))
        #if np.all( a == b):
        if np.array_equal( a, b):
            for i, ia in enumerate( a):
                for j, jb in enumerate( b):
                    if j <= i: continue
                    dist[i,j] = pairwise_dist( ia, jb, dist_fun=dist_fun)
            dist += dist.T
        else:
            for i, ia in enumerate( a):
                for j, jb in enumerate( b):
                    dist[i,j] = pairwise_dist( ia, jb, dist_fun=dist_fun)
        return dist

def scale_data( x, method='range'):
    mmin = x.min (0)
    mmax = x.max (0)
    mmean= x.mean(0)
    mvar = x.var (0)
    if method=='range':
        x_scal = (x - mmin) / (mmax - mmin)
    elif method=='auto':
        x_scal = (x - mmean) / mvar
    elif method=='var':
        x_scal = x / mvar
    elif method=='cent':
        x_scal = x - mmean
    else:
        print( 'Scaling method not recognized, Do not scale')
        x_scal = x.copy()
    return x_scal

class ngas:
    def __init__( self, nn =10, 
                        l_i=30.0, 
                        l_f=0.01, 
                        a_i=0.3, 
                        a_f=0.05, 
                        nk =10,
                        distance_type='eu',
                        print_weights_epochs=False):
        """
        Initialize a Ngas with nn neurons
        INPUT:
        nn                   :: integer - number of neurons
        l_i, l_f, a_i, a_f1  :: float   - initial and final lambda and alpha
        nk                   :: integer - number of lambda neighbors for update
        distance_type        :: str     - name of distance (cb or eu)
        """
        self.num_neurons          = nn
        self._time                = 0
        self._alpha_i             = a_i
        self._alpha_f             = a_f
        self._lamb_i              = l_i
        self._lamb_f              = l_f
        self.nk                   = nk
        self.distance_type        = distance_type
        self.print_weights_epochs = print_weights_epochs
        self.scale_param          = {}

    def load_map( self, weightfile, time=1):
        """
        Read pretrained Ngas from text file
        The scaling method & parameters must be added by the user as
        self.scale_param['method'] = ..
        self.scale_param['min'   ] = ..
        self.scale_param['max'   ] = ..
        self.scale_param['mean'  ] = ..
        self.scale_param['var'   ] = ..
        """
        self.weights = np.genfromtxt( weightfile)
        self.num_neurons, self.n_var = self.weights.shape

    def rescale_inputs( self, scale_method='range', new_params=True):
        """
        Rescale the input data appropriately
        """
        if not hasattr( self, 'scale_method'):
            self.scale_method = scale_method
        if new_params:
            self.scale_param = {}
            mmin = self.train_data.min (0)
            mmax = self.train_data.max (0)
            mmean= self.train_data.mean(0)
            mvar = self.train_data.var (0)
            self.scale_param['method'] = self.scale_method
            self.scale_param['min'   ] = mmin 
            self.scale_param['max'   ] = mmax 
            self.scale_param['mean'  ] = mmean
            self.scale_param['var'   ] = mvar 

        self.original_data = self.train_data.copy()
        if self.scale_method=='range':
            self.train_data = (self.train_data - self.scale_param['min']) /  \
                              (self.scale_param['max'] - self.scale_param['min'])
        elif self.scale_method=='auto':
            self.train_data = (self.train_data - self.scale_param['mean']) /  \
                               self.scale_param['var']
        elif self.scale_method=='var':
            self.train_data = self.train_data / self.scale_param['var']
        elif self.scale_method=='cent':
            self.train_data = self.train_data - self.scale_param['mean']
        else:
            pass

    def scale_back( self, datain=None):
        if datain is None:
            datain = self.weights

        if self.scale_param['method'] == 'range':
            data = ( self.scale_param['max'] - self.scale_param['min'] ) * datain + self.scale_param['min']
        elif self.scale_param['method'] == 'auto':
            data = self.scale_param['var'] * datain + self.scale_param['mean']
        elif self.scale_param['method'] == 'var':
            data = self.scale_param['var'] * datain
        elif self.scale_param['method'] == 'cent':
            data = datain + self.scale_param['mean']
        else:
            data = datain.copy()
        return data

    def train( self
              ,x
              ,scale_method='range'
              ,init_method='data'
              ,upd_method='single'
              ,dist_fun='eucl'
              ,epochs=10
              ,recenter_always=False
              ,endon=None
              ,moved_dist=False):
        """
        Train the map with input data and using the binary distance function
        INPUT:
        np.array(n,m) | x            | input data
        str           | scale_method | choices={'range','cent','var','auto'}
        str           | init_method  | choices={'data','ran'}
        str           | upd_method   | choices={'single',batch}, single is highly recommanded
        lambda fun    | dist_fun     | function to evaluate distance
        int           | epochs       | number of training epochs
        str           | endon        | choices={'traj', 'cent'}
        bool          | moved_dist   | compute how much neurons have moved
        """
        self.train_data            = x
        self.n_samples, self.n_var = self.train_data.shape
        self.scale_method          = scale_method
        self.init_method           = init_method
        self.epochs                = epochs
        self.recenter_always       = recenter_always
        self.endon                 = endon

        # rescale and initialize input data
        self.rescale_inputs()
        self._initialize_net()

        sys.stdout.write("[{0}]".format(" " * epochs))
        sys.stdout.flush()
        sys.stdout.write("\b" * (epochs+1))#back to start of line, after '['

        weights0 = self.weights.copy()

        # train the neurons
        for e in range( epochs):
            if upd_method=='batch':
                self._nd_dist = batch_distance( self.weights, self.train_data, self.distance_type )
                self.calc_deltas()
                self.update_net()
            elif upd_method=='single':
                self.batch_single_upd()
            if self.print_weights_epochs:
                np.savetxt( 'weights{0}.dat'.format(self._time), self.weights)
            self._time += 1

            sys.stdout.write("=")
            sys.stdout.flush()
        sys.stdout.write("]\n")

        if self.endon=='traj':
            # complete training by moving neurons to nearest BMUs (centroids)
            self.complete()        
        elif self.endon=='cent':
            # complete training by moving neurons to centers
            self.recenter_neurons()

        if moved_dist:
            # compute neurons moved distance
            Dw = batch_distance( self.weights, weights0, self.distance_type)
            self.moved_dist = 1.0 / len( weights0) * sum( Dw.min(0))

        # eventually remove redundant neurons
        self.trim_neurons()


    def complete( self):
        """
        Complete training by moving neurons on top of trajectory points:
        move neurons to nearest BMUs
        """
        self.get_BMUs( self.train_data, reverse=True, output=False)
        for i, e in enumerate( self.BMU):
            self.weights[i,:] = self.train_data[e,:]

    def trim_neurons( self):
        """
        Check if there are any redoundant neurons and trim them if it is the 
        case
        """
        mybmus = self.get_BMUs( self.train_data, output=True)
        unique = np.unique( mybmus)
        if unique.shape[0] != self.num_neurons:
            print( 'Neurons redundancy: {0}'.format( self.num_neurons - 
                                                     unique.shape[0]))
            redoundants = np.array( [n for n in np.arange( self.num_neurons) 
                                    if n not in unique])
            print( 'Redoundants to trim:'.format( redoundants))
            self.weights = self.weights[ unique,:]
            self.num_neurons = self.weights.shape[0]

    def batch_single_upd( self): # NEURAL GAS
        """
        Updates the map weights by feeding one data point at a time
        RECOMMENDED UPDATE SCHEME
        """
        data_tmp = self.train_data.copy()
        np.random.shuffle( data_tmp )
        alpha = self._alpha_i * (self._alpha_f / self._alpha_i) ** (float(self._time) / float(self.epochs))
        lamb = self._lamb_i * ( self._lamb_f / self._lamb_i) ** (float(self._time) / float(self.epochs))
        for x in data_tmp:
            D    = batch_distance( np.array( [x]), self.weights, self.distance_type)[0]
            xarg = np.argsort( D)
            xarg = xarg[:self.nk*np.ceil(lamb).astype(int)] if self.nk*lamb < self.num_neurons else xarg
            for i, j in enumerate( xarg):
                self.weights[j] += alpha * np.exp( -i / lamb ) * ( x - self.weights[j] )
        if self.recenter_always:
            self.recenter_neurons()

    def batch_single_upd2( self):
        """
        Alternative single upd
        Less performing apparently
        Neurons are less competitive:
        their neighborhoods do not change within the iteration
        """
        data_tmp = self.train_data.copy()
        np.random.shuffle( data_tmp )
        alpha = self._alpha_i * (self._alpha_f / self._alpha_i) ** (float(self._time) / float(self.epochs))
        lamb = self._lamb_i * ( self._lamb_f / self._lamb_i) ** (float(self._time) / float(self.epochs))
        for neur in self.weights:
            D = batch_distance( np.array( [neur]), data_tmp, self.distance_type)[0]
            xarg = np.argsort( D)
            for i, j in enumerate( xarg):
                xj = data_tmp[j]
                neur += alpha * np.exp( -i / lamb ) * ( xj - neur )
        if self.recenter_always:
            self.recenter_neurons()

    def calc_deltas( self, single_point=None):
        """
        Compute weights updates
        """
        if single_point is not None:
            pass
        self._delta_weights = np.zeros( ( self.num_neurons, self.n_var))
        alpha = ( self.epochs - self._time + 1.0 ) / (self.epochs + 1.0)
        lamb  = ( self.epochs - self._time + 1.0 ) / (self.epochs + 1.0)
        for i, x in enumerate( self.train_data):
            D = self._nd_dist.T[i]
            xarg = np.argsort( D)
            for i, j in enumerate( xarg):
                self._delta_weights[j,:] += alpha * np.exp( -i / lamb ) * ( x - self.weights[j] )

    def recenter_neurons( self):
        """
        Center the neuron to the average of its surrounding samples
        """
        #dmat   = batch_distance( self.weights, self.train_data, self.distance_type ) 
        mybmus = self.get_BMUs( self.train_data, output=True)
        for i, w in enumerate( self.weights):
            neighbors = self.train_data[ np.argwhere( mybmus==i).ravel(),:]
            if neighbors.shape[0] == 0:
                center = w.copy()
            else:
                center = neighbors.mean( axis=0)
            #self.weights[i] = center
            w = center

    def update_net( self):
        """
        Update the network
        optionally center the neurons
        """
        for i, _ in enumerate( self.weights):
            self.weights[i] += self._delta_weights[i] #/ len( self.train_data)
        if self.recenter_always:
            self.recenter_neurons( self)

#    def run_epoch( self):
#        """
#        Run one training epoch
#        """
#        if not hasattr( self, 'train_data'):
#            print( 'You must input training data before start training')
#            return
#        if not hasattr( self, 'weights'):
#            print( 'You must initialize the weights before start training')
#            return
#        if not hasattr( self, 'n_samples'):
#            self.n_samples, self.n_var = self.train_data.shape
#
#        self.epochs = 1
#        self._nd_dist = batch_distance( self.weights, self.train_data, self.distance_type )
#        #self._nn_dist = batch_distance( self.weights, self.weights, self.distance_type )
#        #self.get_BMUs()
#        self.calc_deltas()
#        self.update_net()
#        self._time += 1

    def reset_time( self):
        """
        Reset the _time attribute to 0.
        This increases the learning rate significantly
        """
        self._time = 0
            
    def _initialize_net( self, init_method='data'):
        """
        Initialize the map, either randomly, on random trajectory points,
        of equally spaced trajectory points (default)
        """
        if not hasattr( self, 'init_method'):
            self.init_method = init_method

        if self.init_method == 'ran':
            self.weights = np.random.rand( self.num_neurons, self.n_var)
        elif self.init_method == 'rdata':
            import random
            dataindex = random.sample( np.arange( self.n_samples).tolist(), self.num_neurons)
            #self.weights = self.train_data[ dataindex, :]
            #self.weights = scale_data( self.train_data[ dataindex, :], self.scale_method)
            self.weights = self.train_data[dataindex,:]
        elif self.init_method == 'data':
            chosen = np.round( np.linspace( 0, self.n_samples-1, self.num_neurons)).astype( int)
            #self.weights = scale_data( self.train_data, self.scale_method)[chosen,:]
            self.weights = self.train_data[chosen,:]
        elif self.init_method == 'sphdata':
            self.weights = [self.train_data[0,:],]
            self.sphthr = 3.05
            for i, d in enumerate( self.train_data):
                dist = batch_distance( np.array([d.tolist(),]), 
                                       np.array( self.weights),
                                       self.distance_type).ravel()
                if np.min( dist) > self.sphthr:
                    self.weights.append( d)
            self.weights = np.array( self.weights)
            #self.weights = scale_data( self.weights, self.scale_method)
            self.num_neurons = self.weights.shape[0]
            print( 'Using {0} neurons'.format( self.num_neurons))
        elif self.init_method == 'pca':
            print( 'PCA initialization not implemented yet, using data instead ..')
            self.init_method = 'data'
            self._initialize_net()
        else:
            self.init_method = 'data'
            print( 'Initialization method not recognized. Using "data" method')
            self._initialize_net()


    def get_BMUs( self, xin, reverse=False, output=True, indoutput=True):
        """
        Pass an array of data points in the Ngas coordinates.
        Return Ngas' Best Mathcing Unit for each data point.
        BMU is the closest neuron (Euclidean distance).
        If 'reverse' argument is passed, returns the input point closest to the
        neurons.
        """
        if reverse:
            self.BMU = np.zeros( len( self.weights)).astype( int)
            self._nxin_dist = batch_distance( xin, self.weights, self.distance_type)
        else:
            self.BMU = np.zeros( len( xin))
            self._nxin_dist = batch_distance( self.weights, xin, self.distance_type)
        #for i, di in enumerate( self._nxin_dist.T):
        #    self.BMU[i] = np.argmin( di)
        self.BMU = np.argmin( self._nxin_dist.T, 1)
        if output:
            if indoutput:
                return self.BMU
            else: # full output
                return self.weights[ self.BMU,:] if not reversed else xin[ self.BMU,:]

    def predict( self, x_test):
        """
        Predict a test set of data points
        """
        assert len( x_test.shape) in {1, 2}, 'The input data must a 1D or 2D numpy array'
        self.test_data = x_test
        self.test_pred = self.get_BMUs( self.test_data)
        return self.test_pred
        
if __name__ == '__main__':
    print( 'Nothing to do')
    pass
