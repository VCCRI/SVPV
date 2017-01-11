# if you want to revert to defaults simply re-download this file from github

# pdf viewer, full path to executable or program name
# must work with x11 if using ssh
# default = 'display'
display='display'

# expected library insert size in bp
# not crucial, but can optimise visualisation in some situations
# default = 500
is_len = 500

# expected read length
# not crucial, can optimise visualisation in some situations
# default = 100
read_len = 100

# proportion of sv size to extend window to left and right
# For DEL, DUP, INS and INV only
# default = 1
expansion = 1

# number of bins to split the region into
# will use half this number for single breakppoints
# default = 100
num_bins = 100
