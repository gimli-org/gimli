# -*- coding: utf-8 -*-

import pygimli as g

from math import sqrt, floor, ceil

def unicodeToAscii( text ):
    if type( text ) == unicode:
        return text.encode("iso-8859-1", "ignore")
    else:
        return text

def logDropTol( p, droptol = 1e-3 ):
    '''
    '''
    tmp = g.RVector( p );

    tmp = g.abs( tmp / droptol )
    tmp.setVal( 1.0, g.find( tmp < 1.0 ) )

    #for i, v in enumerate( tmp ):
        #tmp[ i ] = abs( tmp[ i ] / droptol );
        #if tmp[ i ] < 1.0: tmp[ i ] = 1.0;

    tmp = g.log10( tmp );
    tmp *= g.sign( p );
    return tmp;
# def logDropTol

def grange( start, end, dx = 0, n = 0, log = False ):
    '''
        Return RVector from start step-wise filled with dx until end reached [start, end]\n
        Or RVector is filled from start to end with n steps. [start, end] \n
        If nSteps set RVector can optionaly filled with logarithmic spacing.
    '''
    s = float( start )
    e = float( end )
    d = float( dx )
    
    if dx != 0:
        if end < start and dx > 0:
            print "grange: decreasing range but increasing dx, swap dx sign"
            d = -d
        if end > start and dx < 0:
            print "grange: increasing range but decreasing dx, swap dx sign"
            d = -d
        ret = g.asvector( range( int( floor( abs( ( e - s ) / d ) ) + 1 ) ) )
        ret *= d
        ret += s
        return ret;
        
    elif n > 0:
        if not log:
            return grange( start, end, dx = ( e - s ) / ( n - 1 ) )
        else:
            raise Exception( 'not yet implemented.' )
        
    else:
        raise Exception( 'Either dx or nSteps have to be given.' )
    
def diff( v ):
    ''' 
        Return RVector as approximate derivative from v as r[v_1-v_0, v2-v_1,...]
    '''
    r = g.RVector( len( v ) -1 )
    for i in range( len( r ) ):
        r[ i ] = v[ i + 1 ] - v[ i ]
    return r

def xyToLength( x, y ):
    '''
        return RVector of lengths from two RVectors x and y starting from 0 to end
    '''
    ret = g.RVector( len( x ), 0.0 )

    for i in range( len( ret ) -1 ):
        dx = x[ i + 1] - x[ i ]
        dy = y[ i + 1] - y[ i ]

        ret[ i + 1 ] = ret[ i ] + sqrt( dx * dx + dy * dy )
        #ret[ i + 1 ] = ret[ i ] + abs( l[ i + 1 ] - l[ i ] )

    return ret

def getIndex( seq, f ):
    DEPRECATED_SLOW
    idx = [ ]
    if isinstance( seq, g.RVector ):
        for i in range( len ( seq ) ):
            v = seq[ i ]
            if f( v ): idx.append( i )
    else:
        for i, d in enumerate( seq ):
            if f( d ): idx.append( i )
    return idx

def filterIndex( seq, idx ):
    if isinstance( seq, g.RVector ):
        #return seq( idx )
        ret = g.RVector( len( idx ) )
    else:
        ret = range( len( idx ) )

    for i, id in enumerate( idx ):
        ret[ i ] = seq[ id ]
    return ret

def findNearest( x, y, xp, yp, radius = -1 ):
    idx = 0
    minDist = 1e9
    startPointDist = g.RVector( len( x ) )
    for i in range( len( x ) ):
        startPointDist[ i ] = sqrt( ( x[ i ] - xp ) * ( x[ i ] - xp ) + ( y[ i ] - yp ) * ( y[ i ] - yp ) )

        if startPointDist[ i ] < minDist and startPointDist[ i ] > radius:
            minDist = startPointDist[ i ]
            idx = i
    return idx, startPointDist[ idx ]

import itertools

def unique_everseen(iterable, key=None):
    '''List unique elements, preserving order. Remember all elements ever seen.
        http://docs.python.org/library/itertools.html#recipes
    '''
    # unique_everseen('AAAABBBCCDAABBB') --> A B C D
    # unique_everseen('ABBCcAD', str.lower) --> A B C D
    seen = set()
    seen_add = seen.add
    if key is None:
        for element in itertools.ifilterfalse(seen.__contains__, iterable):
            seen_add(element)
            yield element
    else:
        for element in iterable:
            k = key(element)
            if k not in seen:
                seen_add(k)
                yield element

def unique( a ):
    return list( unique_everseen( a ) )

