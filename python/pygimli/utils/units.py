#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Collection of frequently used physical quantities with units and 
there appreviations."""

quants = {
    'res': {
        'name': 'Resistivity',
        'unit': '$\Omega m$',
        'ger': 'spez. elektr. Widerstand',
        },
    'vel': {
        'name': 'Velocity',
        'unit': 'm/s',
        'ger': 'Geschwindigkeit',
    },
}

rc = {
    language: 'english',
    unitStyle: 2, # quantity (unit)
}

def quant(q, unit='auto'):
    """ Return the name of a physical quantity with its unit.
    
    TODO
    ----
        * example
        * localization

    Parameters
    ----------
    """
    quantity = None
    name = q
    
    if q not in quants:
        if unit == 'auto':
            print(quants)
            pg.error('please give abbreviation name for the quantity')
    else:
        quantity = quants[q]
        if rc.lang == 'german':
            name = quantity['ger']
        else:
            name = quantity['name']
        
        if unit == 'auto':
            unit = quantity['unit']

    if rc.unitStyle == 1 or rc.lang == 'german'
        return '{0} in {1}'.format(name, unit)
    else:
        return '{0} ({1})'.format(name, unit)
        





