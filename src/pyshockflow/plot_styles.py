import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm

# Set grid opacity:
grid_opacity = 0.2

# Set font size for axis ticks:
font_axes = 16

# Set font size for axis labels:
font_labels = 24

# Set font size for
font_annotations = 20

# Set font size for plot title:
font_title = 24

# Set font size for plotted text:
font_text = 16

# Set font size for legend entries:
font_legend = 16

# Set font size for colorbar axis label:
font_colorbar = 24

# Set font size for colorbar axis ticks:
font_colorbar_axes = 18

# Set marker size for all line markers:
marker_size_big = 12.5
marker_size = 7.5
marker_size_small = 2
scatter_point_size = 10

# Set the scale for marker size plotted in the legend entries:
marker_scale_legend = 1

# Set line width for all line plots:
heavy_line_width = 3
line_width = 1
medium_line_width = 2
light_line_width = 0.5

# set number of levels in contourf plots
N_levels_coarse = 15
N_levels = 35
N_levels_medium = 70
N_levels_fine = 100
N_fine = 100

# set colormap for contourf plots
color_map = cm.turbo

# number of chars for the banners
total_chars = 100
total_chars_mid = total_chars//2


# Set font size for different elements

# latex preamble
# font_family = 'serif'
# font_name = 'computer Modern'
# plt.rc('text', usetex=True)
# plt.rc('font', family=font_family)
# plt.rc('font', serif=font_name)


plt.rc('font', size=font_text)            # controls default text sizes
plt.rc('axes', titlesize=font_title)       # fontsize of the axes title
plt.rc('axes', labelsize=font_labels)       # fontsize of the x and y labels
plt.rc('xtick', labelsize=font_axes)      # fontsize of the tick labels
plt.rc('ytick', labelsize=font_axes)      # fontsize of the tick labels
plt.rc('legend', fontsize=font_legend)      # legend fontsize
mpl.rcParams['figure.max_open_warning'] = 150
