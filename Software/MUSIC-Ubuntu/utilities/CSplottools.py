#! /usr/bin/env python

from numpy import *


def getPlotElements(idx):
    """
    return line style, marker type, and color for curve. It loops over the
    whole lists of options.
    """
    lineStyleList = ['-', '--', '-.', ':']
    colorList = ['k', 'r', '#009900', 'b', '#990099', '#CC6600', '#009999',
                 '#FF8000', '#00d2ff', '#f12b6e', '#cfc755', '#598254']
    shadowColorList = ['#C0C0C0', '#FFA8A8', '#B3FFB8', '#A4CCFF', '#EF9BFF',
                       '#FFDE80', '#99FFFF', '#FFB266', '#00d2ff', '#fb6eb1',
                       '#cfc755', '#598254']
    MarkerList = ['s', 'o', '^', 'v', 'D', 'p', '*', 'H', '.', ',', '<',
                  '>', '1', '2', '3', '4', 'h', '+', 'x', 'd', '|', '_']
    return (lineStyleList[idx % len(lineStyleList)],
            MarkerList[idx % len(MarkerList)], colorList[idx % len(colorList)],
            shadowColorList[idx % len(shadowColorList)])

