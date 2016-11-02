import plotly.plotly as py
import plotly.graph_objs as go

import numpy as np

py.sign_in('jas_usc', 'sv18e39i18')

layout = {
    'xaxis': {
        'range': [0, 4.5],
        'zeroline': False,
    },
    'yaxis': {        
        'range': [0, 4.5]
    },
    'width': 800,
    'height': 800,
    'shapes': [
        # unfilled circle
        {
            'type': 'circle',
            'xref': 'x',
            'yref': 'y',
            'x0': 1,
            'y0': 1,
            'x1': 3,
            'y1': 3,
            'line': {
                'color': 'rgba(50, 171, 96, 1)',
            },
        },
        # filled circle
        {
            'type': 'circle',
            'xref': 'x',
            'yref': 'y',
            'fillcolor': 'rgba(50, 171, 96, 0.7)',
            'x0': 3,
            'y0': 3,
            'x1': 4,
            'y1': 4,
            'line': {
                'color': 'rgba(50, 171, 96, 1)',
            },
        },
    ]
}



fig = {
    'data': data,
    'layout': layout,
}

py.iplot(fig, filename='shape-circle')














