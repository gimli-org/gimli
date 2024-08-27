#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Stuff for Hankel transformation."""

import numpy as np


def hankelFC(order):
    """Filter coefficients for Hankel transformation.

    10 data points per decade.

    DOCUMENTME .. Author RUB?

    Parameters
    ----------
    order : int

        order=1: NY=+0.5  (SIN)
        order=2: NY=-0.5  (COS)
        order=3: NY=0.0   (J0)
        order=4: NY=1.0   (J1)

    Returns
    -------
    fc : np.array()
        Filter coefficients

    nc0 : int
        fc[nc0] refers to zero argument

    """

    if order == 1:  # sin
        fc = np.array(
            [
                2.59526236e-7,
                3.66544843e-7,
                5.17830795e-7,
                7.31340622e-7,
                1.03322805e-6,
                1.45918500e-6,
                2.06161065e-6,
                2.91137793e-6,
                4.11357863e-6,
                5.80876420e-6,
                8.20798075e-6,
                1.15895083e-5,
                1.63778560e-5,
                2.31228459e-5,
                3.26800649e-5,
                4.61329334e-5,
                6.52101085e-5,
                9.20390575e-5,
                1.30122935e-4,
                1.83620431e-4,
                2.59656626e-4,
                3.66311982e-4,
                5.18141184e-4,
                7.30717340e-4,
                1.03392184e-3,
                1.45742714e-3,
                2.06292302e-3,
                2.90599911e-3,
                4.11471902e-3,
                5.79042763e-3,
                8.20004722e-3,
                1.15192930e-2,
                1.63039133e-2,
                2.28257757e-2,
                3.22249222e-2,
                4.47864328e-2,
                6.27329625e-2,
                8.57059100e-2,
                1.17418314e-1,
                1.53632655e-1,
                1.97717964e-1,
                2.28849849e-1,
                2.40311038e-1,
                1.65409220e-1,
                2.84701476e-3,
                -2.88016057e-1,
                -3.69097406e-1,
                -2.50107514e-2,
                5.71811256e-1,
                -3.92261572e-1,
                7.63280044e-2,
                5.16233994e-2,
                -6.48012082e-2,
                4.89047141e-2,
                -3.26936331e-2,
                2.10539842e-2,
                -1.33862549e-2,
                8.47124695e-3,
                -5.35123972e-3,
                3.37796651e-3,
                -2.13174466e-3,
                1.34513833e-3,
                -8.48749612e-4,
                5.35531006e-4,
                -3.37898780e-4,
                2.13200109e-4,
                -1.34520273e-4,
                8.48765787e-5,
                -5.35535069e-5,
                3.37899801e-5,
                -2.13200365e-5,
                1.34520337e-5,
                -8.48765949e-6,
                5.35535110e-6,
                -3.37899811e-6,
                2.13200368e-6,
                -1.34520338e-6,
                8.48765951e-7,
                -5.35535110e-7,
                3.37899811e-7,
            ],
        )
        nc0 = 40
    elif order == 2:  # cos
        fc = np.array(
            [
                1.63740363e-7,
                1.83719709e-7,
                2.06136904e-7,
                2.31289411e-7,
                2.59510987e-7,
                2.91176117e-7,
                3.26704977e-7,
                3.66569013e-7,
                4.11297197e-7,
                4.61483045e-7,
                5.17792493e-7,
                5.80972733e-7,
                6.51862128e-7,
                7.31401337e-7,
                8.20645798e-7,
                9.20779729e-7,
                1.03313185e-6,
                1.15919300e-6,
                1.30063594e-6,
                1.45933752e-6,
                1.63740363e-6,
                1.83719709e-6,
                2.06136904e-6,
                2.31289411e-6,
                2.59510987e-6,
                2.91176117e-6,
                3.26704977e-6,
                3.66569013e-6,
                4.11297197e-6,
                4.61483045e-6,
                5.17792493e-6,
                5.80972733e-6,
                6.51862128e-6,
                7.31401337e-6,
                8.20645798e-6,
                9.20779729e-6,
                1.03313185e-5,
                1.15919300e-5,
                1.30063594e-5,
                1.45933752e-5,
                1.63740363e-5,
                1.83719709e-5,
                2.06136904e-5,
                2.31289411e-5,
                2.59510987e-5,
                2.91176117e-5,
                3.26704977e-5,
                3.66569013e-5,
                4.11297197e-5,
                4.61483045e-5,
                5.17792493e-5,
                5.80972733e-5,
                6.51862128e-5,
                7.31401337e-5,
                8.20645798e-5,
                9.20779729e-5,
                1.03313185e-4,
                1.15919300e-4,
                1.30063594e-4,
                1.45933752e-4,
                1.63740363e-4,
                1.83719709e-4,
                2.06136904e-4,
                2.31289411e-4,
                2.59510987e-4,
                2.91176117e-4,
                3.26704976e-4,
                3.66569013e-4,
                4.11297197e-4,
                4.61483045e-4,
                5.17792493e-4,
                5.80972733e-4,
                6.51862127e-4,
                7.31401337e-4,
                8.20645797e-4,
                9.20779730e-4,
                1.03313185e-3,
                1.15919300e-3,
                1.30063593e-3,
                1.45933753e-3,
                1.63740362e-3,
                1.83719710e-3,
                2.06136901e-3,
                2.31289411e-3,
                2.59510977e-3,
                2.91176115e-3,
                3.26704948e-3,
                3.66569003e-3,
                4.11297114e-3,
                4.61483003e-3,
                5.17792252e-3,
                5.80972566e-3,
                6.51861416e-3,
                7.31400728e-3,
                8.20643673e-3,
                9.20777603e-3,
                1.03312545e-2,
                1.15918577e-2,
                1.30061650e-2,
                1.45931339e-2,
                1.63734419e-2,
                1.83711757e-2,
                2.06118614e-2,
                2.31263461e-2,
                2.59454421e-2,
                2.91092045e-2,
                3.26529302e-2,
                3.66298115e-2,
                4.10749753e-2,
                4.60613861e-2,
                5.16081994e-2,
                5.78193646e-2,
                6.46507780e-2,
                7.22544422e-2,
                8.03873578e-2,
                8.92661837e-2,
                9.80670729e-2,
                1.07049506e-1,
                1.13757572e-1,
                1.18327217e-1,
                1.13965041e-1,
                1.00497783e-1,
                6.12958082e-2,
                -1.61234222e-4,
                -1.11788551e-1,
                -2.27536948e-1,
                -3.39004453e-1,
                -2.25128800e-1,
                8.98279919e-2,
                5.12510388e-1,
                -1.31991937e-1,
                -3.35136479e-1,
                3.64868100e-1,
                -2.34039961e-1,
                1.32085237e-1,
                -7.56739672e-2,
                4.52296662e-2,
                -2.78297002e-2,
                1.73727753e-2,
                -1.09136894e-2,
                6.87397283e-3,
                -4.33413470e-3,
                2.73388730e-3,
                -1.72477355e-3,
                1.08821012e-3,
                -6.86602007e-4,
                4.33213523e-4,
                -2.73338487e-4,
                1.72464733e-4,
                -1.08817842e-4,
                6.86594042e-5,
                -4.33211523e-5,
                2.73337984e-5,
                -1.72464607e-5,
                1.08817810e-5,
                -6.86593962e-6,
                4.33211503e-6,
                -2.73337979e-6,
                1.72464606e-6,
                -1.08817810e-6,
                6.86593961e-7,
                -4.33211503e-7,
                2.73337979e-7,
                -1.72464606e-7,
            ],
        )
        nc0 = 122
    elif order == 3:  # J0
        fc = np.array(
            [
                2.89878288e-7,
                3.64935144e-7,
                4.59426126e-7,
                5.78383226e-7,
                7.28141338e-7,
                9.16675639e-7,
                1.15402625e-6,
                1.45283298e-6,
                1.82900834e-6,
                2.30258511e-6,
                2.89878286e-6,
                3.64935148e-6,
                4.59426119e-6,
                5.78383236e-6,
                7.28141322e-6,
                9.16675664e-6,
                1.15402621e-5,
                1.45283305e-5,
                1.82900824e-5,
                2.30258527e-5,
                2.89878259e-5,
                3.64935186e-5,
                4.59426051e-5,
                5.78383329e-5,
                7.28141144e-5,
                9.16675882e-5,
                1.15402573e-4,
                1.45283354e-4,
                1.82900694e-4,
                2.30258630e-4,
                2.89877891e-4,
                3.64935362e-4,
                4.59424960e-4,
                5.78383437e-4,
                7.28137738e-4,
                9.16674828e-4,
                1.15401453e-3,
                1.45282561e-3,
                1.82896826e-3,
                2.30254535e-3,
                2.89863979e-3,
                3.64916703e-3,
                4.59373308e-3,
                5.78303238e-3,
                7.27941497e-3,
                9.16340705e-3,
                1.15325691e-2,
                1.45145832e-2,
                1.82601199e-2,
                2.29701042e-2,
                2.88702619e-2,
                3.62691810e-2,
                4.54794031e-2,
                5.69408192e-2,
                7.09873072e-2,
                8.80995426e-2,
                1.08223889e-1,
                1.31250483e-1,
                1.55055715e-1,
                1.76371506e-1,
                1.85627738e-1,
                1.69778044e-1,
                1.03405245e-1,
                -3.02583233e-2,
                -2.27574393e-1,
                -3.62173217e-1,
                -2.05500446e-1,
                3.37394873e-1,
                3.17689897e-1,
                -5.13762160e-1,
                3.09130264e-1,
                -1.26757592e-1,
                4.61967890e-2,
                -1.80968674e-2,
                8.35426050e-3,
                -4.47368304e-3,
                2.61974783e-3,
                -1.60171357e-3,
                9.97717882e-4,
                -6.26275815e-4,
                3.94338818e-4,
                -2.48606354e-4,
                1.56808604e-4,
                -9.89266288e-5,
                6.24152398e-5,
                -3.93805393e-5,
                2.48472358e-5,
                -1.56774945e-5,
                9.89181741e-6,
                -6.24131160e-6,
                3.93800058e-6,
                -2.48471018e-6,
                1.56774609e-6,
                -9.89180896e-7,
                6.24130948e-7,
                -3.93800005e-7,
                2.48471005e-7,
                -1.56774605e-7,
                9.89180888e-8,
                -6.24130946e-8,
            ],
        )
        nc0 = 60
    elif order == 4:  # J1
        fc = np.array(
            [
                1.84909557e-13,
                2.85321327e-13,
                4.64471808e-13,
                7.16694771e-13,
                1.16670043e-12,
                1.80025587e-12,
                2.93061898e-12,
                4.52203829e-12,
                7.36138206e-12,
                1.13588466e-11,
                1.84909557e-11,
                2.85321327e-11,
                4.64471808e-11,
                7.16694771e-11,
                1.16670043e-10,
                1.80025587e-10,
                2.93061898e-10,
                4.52203829e-10,
                7.36138206e-10,
                1.13588466e-9,
                1.84909557e-9,
                2.85321326e-9,
                4.64471806e-9,
                7.16694765e-9,
                1.16670042e-8,
                1.80025583e-8,
                2.93061889e-8,
                4.52203807e-8,
                7.36138149e-8,
                1.13588452e-7,
                1.84909521e-7,
                2.85321237e-7,
                4.64471580e-7,
                7.16694198e-7,
                1.16669899e-6,
                1.80025226e-6,
                2.93060990e-6,
                4.52201549e-6,
                7.36132477e-6,
                1.13587027e-5,
                1.84905942e-5,
                2.85312247e-5,
                4.64449000e-5,
                7.16637480e-5,
                1.16655653e-4,
                1.79989440e-4,
                2.92971106e-4,
                4.51975783e-4,
                7.35565435e-4,
                1.13444615e-3,
                1.84548306e-3,
                2.84414257e-3,
                4.62194743e-3,
                7.10980590e-3,
                1.15236911e-2,
                1.76434485e-2,
                2.84076233e-2,
                4.29770596e-2,
                6.80332569e-2,
                9.97845929e-2,
                1.51070544e-1,
                2.03540581e-1,
                2.71235377e-1,
                2.76073871e-1,
                2.16691977e-1,
                -7.83723737e-2,
                -3.40675627e-1,
                -3.60693673e-1,
                5.13024526e-1,
                -5.94724729e-2,
                -1.95117123e-1,
                1.99235600e-1,
                -1.38521553e-1,
                8.79320859e-2,
                -5.50697146e-2,
                3.45637848e-2,
                -2.17527180e-2,
                1.37100291e-2,
                -8.64656417e-3,
                5.45462758e-3,
                -3.44138864e-3,
                2.17130686e-3,
                -1.36998628e-3,
                8.64398952e-4,
                -5.45397874e-4,
                3.44122545e-4,
                -2.17126585e-4,
                1.36997597e-4,
                -8.64396364e-5,
                5.45397224e-5,
                -3.44122382e-5,
                2.17126544e-5,
                -1.36997587e-5,
                8.64396338e-6,
                -5.45397218e-6,
                3.44122380e-6,
                -2.17126543e-6,
                1.36997587e-6,
                -8.64396337e-7,
                5.45397218e-7,
            ],
        )
        nc0 = 60
    # return (np.reshape(fc, (-1, 1)), nc0)  # (100,) -> (100, 1)
    return fc, nc0
