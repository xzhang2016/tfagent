import seaborn as sns
import matplotlib.pyplot as plt
from datetime import datetime
import os
import pandas as pd
from  matplotlib.colors import LinearSegmentedColormap
import logging
import math

logging.basicConfig(format='%(levelname)s: %(name)s - %(message)s', level=logging.INFO)
logger = logging.getLogger('TFTA-heatmap')


_data_dir = os.path.dirname(os.path.realpath(__file__)) + '/../heatmap/'


def generate_heatmap(filepath):
    """
    Read data from filepath, then generate a heatmap.
    
    parameter
    --------------
    filepath: absolute path to the data
    """
    #sns.set(font_scale=1)
    colors = ['pink', 'blue', 'purple', 'red', 'cyan', 'magenta', 'yellow', 'navy', 'green', 'lime', 'orange']
    try:
        df = pd.read_csv(filepath, sep='\t')
        #check if there's sample type info
        if df.iloc[0][0].lower() == 'sampletype':
            stype = df.iloc[0][1:]
            stype.name = 'Samples'
            df = df.drop(0)
            lut = dict(zip(stype.unique(), colors))
            col_colors = stype.map(lut)
        else:
            col_colors = None
        
        #set index to first column, by default it's gene symbol.
        df = df.set_index('symbol')
        df.index.name = None
        df = df.astype('float')
        
        #calculate zscores for columns/samples
        df = get_zscore(df)
        
        #gene median centered
        #df = median_centered(df, axis=1)
        #calculate zscores for genes/rows. It seems zscore has a better visualization effect than median centered
        df = get_zscore(df, axis=1)
        
        #generate heatmap
        logger.info('Creating the heatmap...')
        if df.shape[0] > 500:
            sns.set(font_scale=0.4)
        g = sns.clustermap(df, method='average', metric='correlation', col_colors=col_colors, \
                           robust=True, cmap=get_cmap())
        #cbar_kws = dict(use_gridspec=False,location="right")
        #whether to show row dendrogram depends on the data size
        if df.shape[0] > 500:
            g.ax_row_dendrogram.set_visible(False)
        
        #change the rotation of ytick label
        plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
        
        #get the reordered row and column index, list
        row_order = g.dendrogram_row.reordered_ind
        col_order = g.dendrogram_col.reordered_ind
        
        #Draw the legend bar for the sidebar
        if col_colors is not None:
            for label in stype.unique():
                g.ax_col_dendrogram.bar(0, 0, color=lut[label], label=label, linewidth=0)
            g.ax_col_dendrogram.legend(title="Samples", loc="best", ncol = math.floor(len(stype.unique())/2))
        
        logger.info('Created heatmap, now try to save to file...')
        filename = save_to_svg(g)
        return filename, row_order, col_order
    except Exception as e:
        logger.info(e)
        return None,None,None

def datetime2str():
    """
    Return a string based on the current date and time.
    """
    now = datetime.now()
    dt_str = now.strftime("%Y%m%dT%H%M%S")
    return dt_str
    
def get_zscore(df, axis=0):
    """
    Calculate z-scores for the rows(axis=1) or the columns(axis=0). 
    Z scores are defined as z = (x - mean)/std.
    
    parameter
    ------------
    df: pandas DataFrame
    """
    mean = df.mean(axis=axis)
    std = df.std(axis=axis)
    if axis:
        #calculate zscore for genes/rows
        df_z = ((df.T - mean)/(std + 0.0001)).T
    else:
        #calculate zscore for samples/columns
        df_z = (df - mean)/(std + 0.0001)
    return df_z
    
def median_centered(df, axis=1):
    """
    Subtract its median from each row if axis=1
    """
    median = df.median(axis=axis)
    if axis:
        df_m = (df.T - median).T
    return df_m
    
def get_cmap():
    #cmap = LinearSegmentedColormap.from_list('gkr',["lime", "black", "orangered"], N=256)
    #or
    #c = ['lime', 'limegreen', 'forestgreen', 'green', 'darkgreen', 'black', 'darkred', 'brown', 'firebrick', 'red', 'orangered']
    c = ['lime', 'limegreen', 'green', 'darkgreen', 'black', 'darkred', 'firebrick', 'red', 'orangered']
    v = [0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0]
    l = list(zip(v,c))
    cmap = LinearSegmentedColormap.from_list('gr',l, N=256)
    return cmap

def save_to_svg(g, data_folder=_data_dir):
    """
    Save g to svg file.
    
    parameter
    ------------
    g: seaborn.matrix.ClusterGrid
    path: str, file path
    """
    if not os.path.isfile(data_folder):
        # Emulate mkdir -p (no error if folder exists)
        try:
            os.mkdir(data_folder)
        except OSError as e:
            if(e.errno != 17):
                logger.error('Could not make the folder: {}.'.format(data_folder))
                return None
    else:
        logger.error('Data path (' + data_folder + ') exists as a file. '
                     'Please rename, remove or change the desired location of the data path.')
        return None
    
    fn = os.path.join(data_folder, datetime2str() + '_heatmap.svg')
    try:
        g.savefig(fn, format='svg')
        #save to pdf only for test
        #fn = os.path.join(data_folder, datetime2str() + '_heatmap.pdf')
        #g.savefig(fn, format='pdf')
        logger.info('Heatmap was saved.')
        return fn
    except Exception:
        logger.info('Could not save file: {}.'.format(fn))
        return None
    
    


