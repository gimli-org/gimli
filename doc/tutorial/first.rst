My first inversion
------------------

Let us start with the very simple inverse problem of fitting a polynomial curve of degree :math:`P`

.. math::

    f(x) = p_0 + p_1 x + \ldots + p_P x^P = \sum\limits_{i=0}^{P} p_i x^i

to some existing data :math:`y`.
The unknown model is the coefficient vector :math:`\m=[p_0,\ldots,p_P]`.
The vectorized function for a vector :math:`\arr{x}=\transpose{[x_1,\ldots,x_N]}` can be written as matrix-vector product

.. _eq:yAx:
.. math::

  \f(\x) = \A \x \quad\mbox{with}\quad \A=\left[ \begin{array}{cccc}
  1 & x_1 & \ldots & x_1^P \\ \vdots & \vdots & \ddots & \vdots \\ 1 & x_N & \ldots & x_N^P
  \end{array} \right] = [ {\bf 1}\quad \x \quad \x^2 \ldots \x^P ] \;.

We set up the modelling operator, i.e. to return :math:`\f(\x)` for given :math:`p_i`, as a class derived from the modelling base class.
The latter holds the main mimic of generating jacobian, gradients by brute force.
The only function to overwrite is \cw{response()}.

*Example file polyfit.py in the directory doc/tutorial/code/polyfit.*

Python is a very flexible language for programming and scripting and has many packages for numerical computing and graphical visualization.
For this reason, we built Python bindings and compiled the library pygimli.
As a main advantage, all classes can be used and derived.
This makes the use of GIMLi very easy for non-programmers.
All existing modelling classes can be used, but it is also easy to create new modelling classes.

We exemplify this by the preceding example.
First, the library must be imported.
To avoid name clashes with other libraries we suggest to import it to an easy name, e.g. by using

>> import pygimli as g

As a result, all gimli objects (classes and functions) can be referred to with a preceding **g.**, e.g.,
g.RVector is the real vector RVector.

Next, the modelling class is derived from ModellingBase, a constructor is defined and the response function is defined.

\begin{lstlisting}[language=python]
import pygimli as g
class FunctionModelling( g.ModellingBase ):
    # constructor
    def __init__( self, nc, xvec, verbose = False  ):
        g.ModellingBase.__init__( self, verbose )
        self.x_ = xvec
        self.nc_ = nc
        self.regionManager().setParameterCount( nc )
    # response function
    def response( self, par ):
        y = g.RVector( self.x_.size(), par[ 0 ] )
        for i in range( 1, self.nc_ + 1 ):
            y += g.pow( self.x_, i ) * par[ i ];
        return y;
    # start model
    def startModel( self ):
        return g.RVector( self.nc_, 0.5 )
\end{lstlisting}

The pygimli library must once be imported (in this case under the name g) and all classes (e.g. modelling operators) can be used by g.classname, e.g. g.RVector is the already known vector of real (double) values.

The main program is very easy then and the code is very similar to C++.
Data are loaded, both forward operator and inversion are created.
Inversion options are set and it the result of run is save to a file.
That's it.

\begin{lstlisting}[language=python]
    xy = g.RMatrix()
    g.loadMatrixCol( xy, datafile );
    # two coefficients and x-vector (first data column)
    f = FunctionModelling( options.np + 1, xy[ 0 ] )
    # initialize inversion with data and forward operator and set options
    inv = g.RInversion( xy[ 1 ], f );
    # constant absolute error of 0.01 (not necessary, only for chi^2)
    inv.setAbsoluteError( 0.01 );
    # the problem is well-posed and does not need regularization
    inv.setLambda( 0 );
    # actual inversion run yielding coefficient model
    coeff = inv.run();
    g.save( coeff, "out.vec" );
\end{lstlisting}

As a main advantage of Python, the actual computations can be easily combined with post-processing or visualization, even building graphical user-interfaces.
In this code example we use matplotlib, a plotting library inside of pylab, a compound of different routines for numerics and plotting, very much comparable to MatLab.

\begin{lstlisting}[language=python]
import pylab as P
P.plot( xy[0], xy[1], 'rx', xy[0], inv.response(), 'b-' )
P.show()
\end{lstlisting}

Similar to C++, command line options can be parsed using the class OptionParser, see the code file.
The output is illustrated for two a synthetic function :math:`y=2.1x+1.1` noisified with Gaussian noise for two different orders in Figure \ref{fig:polyfit}.

\begin{figure}[hbt]%
\includegraphics[width=0.5\columnwidth]{polyfit-n1}\hfill
\includegraphics[width=0.5\columnwidth]{polyfit-n3}%
\caption{Polynomial fit for noisified synthetic data using first order (left) and third order (right) polynomials.}%
\label{fig:polyfit}%
\end{figure}

In the following we continue the description with C++ but all are provided as well in Python without significant code changes.
Main exception is the non-existence of declarations, e.g. the C++ declaration \lstinline|GIMLi::RTransLogLU mytrans( lower, upper );| in Python would be written \lstinline|mytrans = g.RTransLogLU( lower, upper)|.
Whereas it is common to work in name spaces in C++, in python function conflicts easily appear due to lack of declaration.
To avoid that, we recommend importing modules as such, e.g. \lstinline|import pylab as P| and \lstinline|import pygimli as g|.
