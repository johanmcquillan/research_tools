"""This module contains methods to display a command line interface animations."""

import numpy as np
import sys

throbber = {0:  '-',
            1:  '\\',
            2:  '|',
            3:  '/',
            4:  '-',
            5:  '\\',
            6:  '|',
            7:  '/'}
throb_damp = 5E2
bars_total = 20
bar_char = '>'


def update_bar(x, y, bars_done):
    fraction = float(x) / y
    bars = fraction * bars_total
    if bars >= bars_done:
        bars_done = int(np.round(bars))
        sys.stdout.write('\r  [{:<{}}] '.format(bar_char*bars_done, bars_total))
        sys.stdout.write(' {:>3.0f}%'.format(np.ceil(fraction * 100)))
        sys.stdout.flush()
    return bars_done
