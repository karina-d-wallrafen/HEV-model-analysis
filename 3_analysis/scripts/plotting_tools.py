import matplotlib.pyplot as plt
import matplotlib._mathtext as mathtext

def set_fontsize(base_fontsize=16):
    fontsize = base_fontsize
    plt.rcParams.update({
        'font.size': fontsize,
        'axes.titlesize': fontsize * 1.125,
        'axes.labelsize': fontsize,
        'xtick.labelsize': fontsize,
        'ytick.labelsize': fontsize,
        'legend.fontsize': fontsize,
        'font.family': "Arial",
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'Arial',
        'mathtext.it': 'Arial:italic',
        'mathtext.bf': 'Arial:bold',
    })
    mathtext.FontConstantsBase.sup1 = 0.4

plt.style.use('default')

set_fontsize()

teal = "#21918C"
lightblue = "#528DBE"
darkblue = "#1F4E79"
lightgreen = "#4DAF4A"
darkgreen = "#276419"
purple = "#5E2A70"
orange = "#E66101"
maroon = "#7F0000"
lightest_grey = '#eeeeee'
lightgrey = '#cccccc'
medium_grey = '#a2a2a2'
darkgrey = '#595959'
darkest_grey = '#3d3d3d'

