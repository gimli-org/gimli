My first inversion
------------------

.. literate:: first.py

Similar to C++, command line options can be parsed using the class OptionParser, see the code file.
The output is illustrated for two a synthetic function :math:`y=2.1x+1.1` noisified with Gaussian noise for two different orders in Figure \ref{fig:polyfit}.



In the following we continue the description with C++ but all are provided as well in Python without significant code changes.
Main exception is the non-existence of declarations, e.g. the C++ declaration \lstinline|GIMLi::RTransLogLU mytrans( lower, upper );| in Python would be written ref::`g.RTransLogLU`

Das ist Beispielcode::    

    mytrans = g.RTransLogLU( lower, upper)

::    

    mytrans = g.RTransLogLU( 'lower', upper)

Whereas it is common to work in name spaces in C++, in python function conflicts easily appear due to lack of declaration.
To avoid that, we recommend importing modules as such, e.g. \lstinline|import pylab as P| and \lstinline|import pygimli as g|.
