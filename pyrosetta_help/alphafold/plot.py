__all__ = ['make_pae_plot']

import numpy as np

def make_pae_plot(errors: np.ndarray) -> 'plotly.import graph_objs.Figure':
    """
    Make AlphaFold2-EBI like PAE plot
    """
    import plotly.express as px
    fig = px.imshow(errors, color_continuous_scale=[(0, 'green'), (1, 'white')])
    return fig
