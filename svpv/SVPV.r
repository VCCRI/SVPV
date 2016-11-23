# Container for plot parameters
PlotParams <- function(args) {
  this <- list(
    depth = ('-d' %in% args),
    orphaned = ('-or' %in% args),
    inverted = ('-v' %in% args),
    samestrand = ('-ss' %in% args),
    secondary = ('-se' %in% args),
    supplementary = ('-su' %in% args),
    hardclipped = ('-hc' %in% args),
    ins = ('-i' %in% args),
    refgene = ('-r' %in% args),
    svAF = ('-af' %in% args),
    legend = ('-l' %in% args)
  )
  class(this) <- append(class(this), "PlotParams")
  return (this)
}

# Container for data for a given sample
Sample <- function(folder, sample, label = FALSE) {
  depths <- read.delim(paste0(folder, sample, '/depths.tsv'), header = TRUE,  sep = '\t')
  bin_size <- depths$bin[2] -  depths$bin[1]
  num_bins <- nrow(depths)
  xlims <- c(depths$bin[1], (depths$bin[nrow(depths)] + bin_size))
  #read in SV calls for this sample
  svs <- read.delim(paste0(folder, sample, '/svs.tsv'), header = TRUE, sep = '\t')
  #convert to list, separating calls from different vcfs
  svs <- split.data.frame(svs, svs$vcf)
  svTracks <- lapply(svs, function(x) get_tracks(x$start, x$end))
  #read in foward and reverse inserts
  fwd_ins <- read.table(paste0(folder, sample, '/fwd_ins.tsv'), fill=TRUE, sep="\t", stringsAsFactors = FALSE, header = FALSE)
  rvs_ins <- read.delim(paste0(folder, sample, '/rvs_ins.tsv'), fill=TRUE, sep="\t", stringsAsFactors = FALSE, header = FALSE)
  #convert inserts from strings to numerics
  fwd_ins <- cbind(fwd_ins[,1], lapply(fwd_ins[,2], function(x) as.numeric(unlist(strsplit(x, ',')))))
  rvs_ins <- cbind(rvs_ins[,1], lapply(rvs_ins[,2], function(x) as.numeric(unlist(strsplit(x, ',')))))
  #read aln_stats
  aln_stats <- read.delim(paste0(folder, sample, '/aln_stats.tsv'), header = TRUE,  sep = '\t')
  aln_stats_bin_size <- aln_stats$bin[2] - aln_stats$bin[1]
  split <- FALSE
  for (i in 3:nrow(aln_stats)){
    if (aln_stats_bin_size != aln_stats$bin[i] - aln_stats$bin[i-1]){
      split = i
    }
  }
  ins_ylim <- max(max(unlist(sapply(fwd_ins, function(x)  estimate_upper_bound(x)))),  max(unlist(sapply(rvs_ins, function(x) estimate_upper_bound(x)))))
  
  this <- list(
    Name = sample,
    Depths = depths,
    Bin_size = bin_size,
    Num_bins = num_bins,
    Xlims = xlims,
    SVs = svs,
    SVtracks = svTracks,
    Fwd_ins = fwd_ins,
    Rvs_ins = rvs_ins,
    Aln_stats = aln_stats,
    Split = split,
    Ins_ylim = ins_ylim
  )
  class(this) <- append(class(this), "Sample")
  return(this)
}
#container for annotations
Annotations <- function(folder) {
  #check if refgenes file exists, if so read
  genes_file <- paste0(folder, 'refgene.tsv')
  if (file.exists(genes_file)) {
    genes <- read.delim(genes_file, header = TRUE, sep = '\t')
  } else {
    genes = NULL
  }
  #check if SV_AF file exists, if so read
  SV_AF_file = paste0(folder, 'SV_AF.tsv')
  if (file.exists(SV_AF_file)) {
    SV_AF <- read.delim(SV_AF_file, header = TRUE, sep = '\t')
    SV_AF <- split.data.frame(SV_AF, SV_AF$vcf)
    AF_tracks <-
      lapply(SV_AF, function(x)
        get_tracks(x$start, x$end))
  } else {
    SV_AF = NULL
    AF_tracks = NULL
  }
  
  this <- list(
    Genes = genes,
    SV_AF = SV_AF,
    AF_tracks = AF_tracks
  )
  class(this) <- append(class(this), "Annotations")
  return(this)
}

# add a border around a given plot area
add_border <- function(xlims, ylims, lwd = 0.5) {
  rect(xlims[1], ylims[1], xlims[2], ylims[2], lwd = lwd)
}

# horizontal line separator between samples
separator <- function() {
  empty_plot(c(0, 1))
  par(xpd = NA)
  abline(h = 0.5, lwd = 4, col = 'gray25')
  par(xpd = FALSE)
}

# create an empty plot
empty_plot <- function(xlim,ylim = c(0, 1),type = 'n',bty = 'n', xaxt = 'n', yaxt = 'n', ylab = '', xlab = '') {
    plot(1, type = type, ylim = ylim, xlim = xlim, bty = bty, xaxt = xaxt, yaxt = yaxt, ylab = ylab, xlab = xlab)
}

# add axis to a plot
add_position_axis <- function(xlims, side) {
  units <- get_units(xlims[2] - xlims[1])
  empty_plot(xlims / units$val)
  mtext( paste0('Position (', units$sym, ')'), side = side, line = -1, cex = 0.85
  )
  axis(side = side, line = -3)
}

# plot read depth and mapping quality
plot_depth <- function(depth, xlims) {
  par(las = 1)
  bin_size = depth$bin[2] - depth$bin[1]
  ylims = c(0, 1.2 * max(depth$total - depth$mapQ0 - depth$mapQltT, na.rm = TRUE))
  empty_plot(xlims, ylim = ylims)
  title(ylab = 'Depth\n(reads/ bp)', line=2)
  add_border(xlims, ylims)
  # add depth$total depth
  rect(depth$bin, c(0), (depth$bin + bin_size), depth$total, col = '#74C476')
  # add depth$mapQ0 to depth
  rect(depth$bin, depth$total, (depth$bin + bin_size), (depth$total - depth$mapQ0), col = 'white')
  # add depth$mapQltT to depth
  rect(depth$bin, (depth$total - depth$mapQ0), (depth$bin + bin_size), depth$total - (depth$mapQ0 + depth$mapQltT), col = 'khaki1')
  axis(2, tick = TRUE, labels = TRUE, line = -1)
}

# add the legend
add_legend <- function() {
  # create a plot with room for four legends: depth, inserts, mapping stats, svtype/freq
  empty_plot(c(0, 4), ylim = c(0, 1))
  add_border(c(0, 1), c(0, 1))
  add_border(c(1, 2), c(0, 1))
  add_border(c(2, 3), c(0, 1))
  add_border(c(3, 4), c(0, 1))
  par(font = 2)
  text(c(0.5, 1.5, 2.5, 3.5), 1,  pos = 1, labels = c("Read MapQ", "Inferred Insert Size", "Mapping Stats", "SV Allele Frequency"))
  par(font = 1)
  # constants
  bottom = 0.20
  top = bottom + 0.20
  # depth legend
  rect(c((1 / 6 - 0.1), (3 / 6 - 0.1), (5 / 6 - 0.1)), bottom, c((1 / 6 + 0.1), (3 / 6 + 0.1), (5 / 6 + 0.1)), top,  col = c('#74C476', 'khaki1', 'gray95'))
  text( c((1 / 6), (3 / 6), (5 / 6)), c(top + 0.1), pos = 3, labels = c(">= 30", "< 30", "= 0") )
  
  # inferred insert size legend
  
  text(1.5, top + 0.25, labels = c("Proportion in position x\nwith mapping distance y"))
  rect(1.15 + 0.7 / 10 * (0:9), bottom, 1.15 + 0.7 / 10 * (1:10), top, col=insert_size_pallete(10)[(1:10)])
  text(c(1.18, 1.82), bottom - 0.08, as.character(c(0, 1)))
  
  # Mapping stats legend 
  text(2.5, top + 0.1, pos = 3, labels = c("proportion of reads in position x"))
  rect(2.15 + 0.7 / 10 * (0:9), bottom, 2.15 + 0.7 / 10 * (1:10), top, col=aln_stats_pallete(10)[(1:10)])
  text(c(2.18, 2.82), bottom - 0.08, as.character(c(0, 1)))
  
  # SV AF legend
  sv_types <- c("DEL", "DUP", "CNV", "INV")
  height = 0.13
  top = bottom + 0.15 * length(sv_types)
  par(family='mono', font = 2)
  for (i in 1:length(sv_types)){
    rect(3.20 + 0.7 / 10 * (0:9), bottom + (i-1) * height, 3.20 + 0.7 / 10 * (1:10), bottom + (i) * height, col=sapply(1:10, function(x) get_sv_col(sv_types[i], x/10)))
    text(3.20, bottom + (i-0.5) * height, sv_types[i], pos=2)
  }
  par(family='sans', font = 1)
  text(c(3.23, 3.87), bottom - 0.08, as.character(c(0, 1)))
}

aln_stats_pallete <- function(n){
  return(colorRampPalette(c("gray95", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603", "#7F2704"))(n))
}
# plot alignment stats
plot_aln_stats <- function(total, numerator, label, split, spacer=2) {
    if (split) {end = length(total)+4*spacer} else { end = length(total)}
    empty_plot(c(0, end))
    par(las = 1)
    mtext(label, side = 2, line = -1,  cex = 0.75)
    #brewer YlGnBu pallete
    if (split){
      rect(spacer:(spacer+split-2), c(0), (spacer+1):(spacer+split-1), c(1), col=aln_stats_pallete(20)[(19 * numerator[1:(split-1)] / total[1:(split-1)]) + 1], border =NA)
      add_border(c(spacer, spacer+split-1), c(0, 1))
      rect((3*spacer+split-1):(3*spacer+length(total)-1), c(0), (3*spacer+split):(3*spacer+length(total)), c(1), col=aln_stats_pallete(20)[(19 * numerator[split:length(total)] / total[split:length(total)]) + 1], border =NA)
      add_border(c(3*spacer+split-1, 3*spacer+length(total)), c(0, 1))
    } else {
      rect(0:(end-1), c(0), 1:end, c(1), col=aln_stats_pallete(20)[(19 * numerator[1:end] / total[1:end]) + 1], border =NA)
      add_border(c(0, end), c(0, 1))
    }
}

# returns a list of the form (val=, sym=)
get_units <- function(num_bp) {
  if (num_bp < 1000) {
    return(list(val = 1, sym = 'bp'))
  } else if (num_bp < 1000000) {
    return(list(val = 1000, sym = 'kbp'))
  } else if (num_bp < 1000000000) {
    return(list(val = 1000000, sym = 'Mbp'))
  } else {
    return(list(val = 1000000000, sym = 'Gbp'))
  }
}

# return the lowest insert size in the highest top 15%
# simple heuristic for avoiding outliers (assumes outliers are less than 15% abundant, and true inserts are at least 15% abundant)
estimate_upper_bound <- function(ins) {
  return(1.1 * sort(ins)[floor(length(ins) * 0.85)])
}

insert_size_pallete <- function(n){
  return(colorRampPalette(c("gray95", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B"))(n))
}

# plots the inserts in specifiend interval
plot_binned_inserts <- function(binned_inserts, num_y_bins, split, spacer=2){
  if (split) {end = nrow(binned_inserts) + 4*spacer} else {end = nrow(binned_inserts)}
  empty_plot(c(0,end), ylab = '' ,ylim = c(0, num_y_bins + 1))
  if (split){
    for (i in 1:(num_y_bins + 1)) {
      rect(spacer:(spacer+split-2), i - 1, (spacer+1):(spacer+split-1), i, col=insert_size_pallete(25)[24 * binned_inserts[1:(split-1), i] + 1],  border = NA)
      rect((3*spacer+split-1):(3*spacer+nrow(binned_inserts)-1), i - 1, (3*spacer+split):(3*spacer+nrow(binned_inserts)), i, col=insert_size_pallete(25)[24 * binned_inserts[(split:nrow(binned_inserts)), i] + 1],  border = NA)
      add_border(c(spacer, spacer+split-1), c(0, num_y_bins))
      add_border(c(spacer, spacer+split-1), c(num_y_bins, num_y_bins + 1))
      add_border(c(3*spacer+split-1, 3*spacer+nrow(binned_inserts)), c(0, num_y_bins))
      add_border(c(3*spacer+split-1, 3*spacer+nrow(binned_inserts)), c(num_y_bins, num_y_bins + 1))
    }
  } else {
    for (i in 1:(num_y_bins + 1)) {
      rect((1:nrow(binned_inserts))-1, i - 1, (1:nrow(binned_inserts)), i, col=insert_size_pallete(25)[24 * binned_inserts[(1:nrow(binned_inserts)), i] + 1],  border = NA)
      add_border(c(0, end), c(0, num_y_bins))
      add_border(c(0, end), c(num_y_bins, num_y_bins + 1))
    }
  }
}

plot_insert_sizes <- function(fwd_ins, rvs_ins, ylim, split, num_y_bins = 10) {
  # divide into 10 bins spaced equally between 0 and ylim
  ybin_size <- ylim / num_y_bins
  # create an extra bin to store anythin larger than ylim
  fwd_bins <- matrix(nrow = nrow(fwd_ins), ncol = (num_y_bins + 1))
  rvs_bins <- matrix(nrow = nrow(rvs_ins), ncol = (num_y_bins + 1))
  # get counts for each bin
  for (i in 1:num_y_bins) {
    fwd_bins[, i] = sapply(fwd_ins[,2], function(x) sum((((i - 1) * ybin_size)  <= x) & (x < ((i) * ybin_size))))
    rvs_bins[, i] = sapply(rvs_ins[,2], function(x) sum((((i - 1) * ybin_size)  <= x) & (x < ((i) * ybin_size))))
  }
  fwd_bins[, (num_y_bins + 1)] = sapply(fwd_ins[,2], function(x) sum(ylim <= x))
  rvs_bins[, (num_y_bins + 1)] = sapply(rvs_ins[,2], function(x) sum(ylim <= x))
  
  # convert to proportions
  fwd_bins = fwd_bins / rowSums(fwd_bins)
  # fwd_bins[is.nan(fwd_bins)] <- 0
  rvs_bins = rvs_bins / rowSums(rvs_bins)
  # rvs_bins[is.nan(rvs_bins)] <- 0
  
  # organise sensible units for ticks on plot
  units <- get_units(ylim / 2)
  ylim <- ylim / units$val
  mid <- round(ylim / 2)
  interval <- round(mid * 2 / 3)
  ticks_at <- c((mid - interval), mid, (mid + interval))

  # plot binned foward inserts
  par(las = 1)
  plot_binned_inserts(fwd_bins, num_y_bins, split)
  axis(2, at = 10 * ticks_at / ylim,  labels = as.character(ticks_at),  line = -1)
  title(ylab=paste0('forward\nmapping\ndistance\n (', units$sym, ')'), line=1.5)
  par(xpd = NA)
  text(0, num_y_bins + 0.5, labels =">", cex = 0.85, pos = 2)
  
  # plot binned reverse inserts
  plot_binned_inserts(rvs_bins, num_y_bins, split)
  axis( 2, at = 10 * ticks_at / ylim,  labels = as.character(ticks_at), line = -1)
  title(ylab=paste0('reverse\nmapping\ndistance\n (', units$sym, ')'), line=1.5)
  text(0, num_y_bins + 0.5, labels =">", cex = 0.85, pos = 2)
  par(xpd = FALSE)
}
# plot structural variants
plot_svs <- function(svs, xlims, tracks, AF=TRUE) {
  empty_plot(xlims)
  add_border(xlims,c(0,1))
  mtext(svs$vcf, side = 2, line = -1, cex = 0.8)
  if (is.null(svs)) {
    text(0.5 * (xlims[1] + xlims[2]), 0.5, labels = "None")
  }
  # get mapping of svs to tracks to ensure no overlap
  scale = 1 / max(tracks)
  par(las = 1)
  par(font = 2)
  for (i in 1:nrow(svs)) {
    # create a rectange covering each sv call
    start <- max(xlims[1], svs$start[i])
    end <- min(xlims[2], svs$end[i])
    x_prop <- (end - start) / (xlims[2] - xlims[1])
    bottom <- ((tracks[i] - 1) * scale)
    top <- ((tracks[i]) * scale)
    spacer = 0.1*(top-bottom)
    # label the sv, if it is of sufficient length to not overlap bounds
    len <- svs$end[i] - svs$start[i] + 1
    units <- get_units(len)

    if (!AF){
      rect(start, bottom + spacer, end, top - spacer, col=get_sv_col(svs$svtype[i], 0.8*(0.5*grepl('1', svs$gt[i]) + 0.5*grepl('1/1', svs$gt[i]))), border=get_sv_col(svs$svtype[i], 1), lwd=2)
      if (x_prop > 1/5){
        text(  0.5 * (max(xlims[1], svs$start[i]) + min(xlims[2], svs$end[i])), ((tracks[i] - 0.5) * scale), labels = paste(svs$svtype[i], ':', svs$gt[i], ':', as.character(round(len/units$val, digits=2)), " ", units$sym))
      } else {
        text(  0.5 * (max(xlims[1], svs$start[i]) + min(xlims[2], svs$end[i])), ((tracks[i] - 0.5) * scale), labels = svs$gt[i])
      }
    } else {
      rect(start, bottom + spacer, end, top - spacer, col=get_sv_col(svs$svtype[i], as.numeric(svs$MAF[i])), border=get_sv_col(svs$svtype[i], 1), lwd=2)
      if (x_prop > 1/5){
        text(0.5 * (start + end), 0.5 * (top + bottom), labels = paste0(svs$svtype[i], ' : AF = ', as.character(round(as.numeric(svs$MAF[i]), digits = 3)), ' : ', as.character(round(len/units$val, digits=2)), " ", units$sym))
      } else if (x_prop > 1/10){
        text(0.5 * (start + end), 0.5 * (top + bottom), labels = paste0('AF = ', as.character(round(as.numeric(svs$MAF[i]), digits = 3))))
      }
    }
  }
  par(font = 1)
}

# for getting sv col based on genotype
gt_to_intensity <- function(gt){
  if (grepl("0/1", gt)){
    return(0.3)
  } else if (grepl("0/1", gt)){
    return(0.8)
  } else {
    return(0)
  }
}
# return a colour for a given SV type
# intensity is a value between 0 and 1, if intensity is zero white is alwaya returned
get_sv_col <- function(type, intensity) {
  if ((is.na(intensity)) | (is.nan(intensity))){
    return('white')
  } else if (intensity == 0){
    return('white')
  } else if ((grepl('DEL', type))) {
      return(colorRampPalette(c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15"))(100)[20 + floor(80 * intensity)])
  } else if ((grepl('DUP', type))) {
      return(colorRampPalette(c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C"))(100)[20 + floor(80 * intensity)])
  } else if ((grepl('INV', type))) {
      return(colorRampPalette(c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C"))(100)[20 + floor(80 * intensity)])
  } else if ((grepl('CNV', type))) {
    return(colorRampPalette(c("#FCFBFD", "#EFEDF5", "#DADAEB", "#BCBDDC", "#9E9AC8", "#807DBA", "#6A51A3", "#54278F"))(100)[20 + floor(80 * intensity)])
  } else {
    return(colorRampPalette(c("#FFFFFF", "#F0F0F0", "#D9D9D9", "#BDBDBD", "#969696", "#737373", "#525252", "#252525"))(100)[20 + floor(80 * intensity)])
  }
}

# plot specified tracks for a given sample
plot_sample <- function(sample, plot_params, ins_ylim) {
  #plot sample title
  empty_plot(sample$Xlims)
  par(font = 2)
  text(sample$Xlims[1], 0.5, paste("Sample:", sample$Name), cex = 1, pos = 4)
  par(font = 1)
  #plot sample SVs
  par(las=1)
  for (i in 1:length(sample$SVs)) {
    plot_svs(sample$SVs[[i]], sample$Xlims, sample$SVtracks[[i]], AF=FALSE)
  }
  # plot sample depth
  if (plot_params$depth) { plot_depth(sample$Depths, sample$Xlims)}
  # plot zoom details
  if (sample$Split){
    add_zoom_detail(sample$Xlims, sample$Aln_stats$bin[1], sample$Aln_stats$bin[sample$Split-1], sample$Aln_stats$bin[sample$Split], sample$Aln_stats$bin[length(sample$Aln_stats$bin)], length(sample$Aln_stats$bin))
  }
  # plot sample insert sizes
  if (plot_params$ins) {
    plot_insert_sizes(sample$Fwd_ins, sample$Rvs_ins, ins_ylim, sample$Split)
  }
  # plot remaining tracks
  if (plot_params$hardclipped) {
    plot_aln_stats( sample$Aln_stats$reads, sample$Aln_stats$hardclipped, 'hardclipped', sample$Split)
  }
  if (plot_params$secondary) {
    plot_aln_stats( sample$Aln_stats$reads, sample$Aln_stats$secondary, 'secondary', sample$Split)
  }
  if (plot_params$supplementary) {
    plot_aln_stats( sample$Aln_stats$reads, sample$Aln_stats$supplementary, 'supplementary', sample$Split)
  }
  if (plot_params$orphaned) {
    plot_aln_stats( sample$Aln_stats$reads, sample$Aln_stats$orphaned, 'orphaned', sample$Split)
  }
  if (plot_params$inverted) {
    plot_aln_stats( sample$Aln_stats$reads, sample$Aln_stats$inverted, 'inverted', sample$Split)
  }
  if (plot_params$samestrand) {
    plot_aln_stats( sample$Aln_stats$reads, sample$Aln_stats$samestrand, 'samestrand', sample$Split)
  }
  # plot zoom details
  if (sample$Split){
    add_zoom_detail(sample$Xlims, sample$Aln_stats$bin[1], sample$Aln_stats$bin[sample$Split-1], sample$Aln_stats$bin[sample$Split], sample$Aln_stats$bin[length(sample$Aln_stats$bin)], length(sample$Aln_stats$bin), axes=TRUE)
  }
  # add separator
  separator()
}

plot_details <- function(bin_size, num_bins) {
  mtext( paste("Bin size: ", as.character(bin_size), "    Num bins: ", as.character(num_bins), "    Date: ", as.character(Sys.Date()) ), side = 1, line = 0, adj = 0, cex = 0.65 )
}
# graphical representation of linear transformation between depth plot and zoomed in inserts size and aln stats plots
add_zoom_detail <- function(xlims, start_1, end_1, start_2, end_2, num_bins, col='black', spacer=2, axes = FALSE){
  range = xlims[2] - xlims[1]
  zoom_start_1 = xlims[1] + (spacer/(num_bins+4*spacer))*range
  zoom_end_1 = xlims[1] + ((spacer + 0.5*num_bins)/(num_bins+4*spacer))*range
  zoom_start_2 = xlims[1] + ((3*spacer+0.5*num_bins)/(num_bins+4*spacer))*range
  zoom_end_2 = xlims[1] + ((3*spacer+num_bins)/(num_bins+4*spacer))*range
  empty_plot(xlims)
  if (axes){
    # add axes for the zoomed regions
    units <- get_units(end_1-start_1)
    at = (zoom_start_1 + (zoom_end_1-zoom_start_1)*c((1/6),(1/2),(5/6)))
    labels = ((start_1 + (end_1-start_1)*c((1/6),(1/2),(5/6)))/units$val)
    axis(side=1, at=at, labels=paste(as.character(round(labels, digits=1)), units$sym), line=-2)
    at = (zoom_start_2 + (zoom_end_2-zoom_start_2)*c((1/6),(1/2),(5/6)))
    labels = ((start_2 + (end_2-start_2)*c((1/6),(1/2),(5/6)))/units$val)
    axis(side=1, at=at, labels=paste(as.character(round(labels, digits=1)), units$sym),line=-2)
  } else {
    # show the level of zoom
    segments(c(start_1, end_1, start_2, end_2), c(0.7), c(start_1, end_1, start_2, end_2), c(2), col=col, lwd=2)
    segments(c(start_1, end_1, start_2, end_2), c(0.7), c(zoom_start_1, zoom_end_1, zoom_start_2, zoom_end_2), c(0.3), col=col, lwd=2)
    segments(c(zoom_start_1, zoom_end_1, zoom_start_2, zoom_end_2), c(0.3), c(zoom_start_1, zoom_end_1, zoom_start_2, zoom_end_2), c(-1), col=col, lwd=2)
  }
}

#returns an assignment to tracks for a set of regions such that there are no overlaps
get_tracks <- function(starts, ends) {
  # assume that starts and ends are of same length
  # assume also that they are sorted by lowest start first
  tracks <- vector("integer", length = length(starts))
  tracks[1] = 1
  if (length(starts) >= 2) {
    for (i in 2:length(starts)) {
      for (j in 1:i) {
        overlap = FALSE
        if (j %in% tracks) {
          check = which(tracks %in% j)
          for (k in 1:length(check)) {
            if (starts[i] <= ends[check[k]]) {
              overlap = TRUE
              break
            }
          }
        }
        if (!overlap) {
          tracks[i] = j
          break
        }
      }
    }
  }
  return(tracks)
}

plot_refgenes <- function(refgenes, xlims) {
  empty_plot(xlims)
  plot_range <- xlims[2] - xlims[1]
  # if no refgene annotation in region don't plot it
  if (is.null(refgenes)) {
    text(0.5 * (xlims[1] + xlims[2]), 0.5, labels = "None")
  } else {
    # ensure no genes are plotted overlapping by assigning those that do overlap to separate tracks in a greedy fashion
    # since refgenes should already be sorted by start position this is fairly straightforward
    tracks <- 1:nrow(refgenes)
    scale = 1/max(tracks)
    for (i in 1:(nrow(refgenes))) {
      # plot thin rectangle for whole length of transcript
      plot_start = max(xlims[1], refgenes$txStart[i])
      plot_end = min(xlims[2], refgenes$txEnd[i])
      plot_len = plot_end - plot_start + 1
      total_len = refgenes$txEnd[i] - refgenes$txStart[i] + 1
      fwd = ("+" == as.character(refgenes$strand[i]))
      if (fwd) {dir = '->'} else {dir = '<-'}
      segments(xlims[1], (tracks[i]- 0.5) * scale, xlims[2], (tracks[i]- 0.5) * scale, col='gray50')
      rect(plot_start, ((tracks[i] - 0.8) * scale), plot_end, ((tracks[i]-0.2) * scale), col = '#74C476', border='gray50')
      units <- get_units(plot_len)
      par(font=2)
      #l abel the gene
      par(xpd=NA)
      text(xlims[1], ((tracks[i]- 0.5) * scale), labels=paste(refgenes$name2[i], dir), pos=2)
      text(xlims[2], ((tracks[i]- 0.5) * scale), labels=paste(as.character(round(100*plot_len/total_len, digits=0)), '%'), pos=4)
      par(xpd=FALSE)
      # get the exon starts and ends
      starts = as.numeric(strsplit(as.character(refgenes$exonStarts[i]), ',')[[1]])
      ends = as.numeric(strsplit(as.character(refgenes$exonEnds[i]), ',')[[1]])
      fwd = ("+" == as.character(refgenes$strand[i]))
      min_exon_label_dist = 0.02 * (xlims[2] - xlims[1])
      last_exon_labelled = NA
      # plot exons
      for (j in 1:length(starts)) {
        # check if exon is within plot limits
        if ((ends[j] < xlims[1]) | (starts[j] > xlims[2])) { next }
        rect(max(xlims[1], starts[j]), ((tracks[i] - 0.925) * scale),  min(xlims[2], ends[j]), ((tracks[i] - 0.075) * scale), col = '#6BAED6', border='gray50')
      }
      for (j in 1:length(starts)) {
        if ((ends[j] < xlims[1]) | (starts[j] > xlims[2])) { next }
        if (fwd) {num = j} else {num = length(starts) - j + 1}
        label_pos = 0.5 * (max(xlims[1], starts[j]) + min(xlims[2], ends[j]))
        # ensure enought distance between exon labels and gene name label before annotating
        if (is.na(last_exon_labelled) | (label_pos - last_exon_labelled) > min_exon_label_dist){
          text(0.5 * (max(xlims[1], starts[j]) + min(xlims[2], ends[j])), ((tracks[i]-0.5) * scale), labels = as.character(num))
          last_exon_labelled = label_pos
        }
      }
    }
    par(font=1)
    }
  }

get_plot_layout <- function(plot_params, annotations, num_samples, vcfs_per_sample, split, max=200) {
    # note: order of plots to be as implied here
    # top x-axis
    heights <- c(3)
    # title
    heights <- c(heights,1)
    # separator
    heights <- c(heights, 1)
    for (i in 1:num_samples){
      # sample title
      heights <- c(heights, 1.5)
      # add in vcf plots
      for (j in 1:length(vcfs_per_sample[[i]])){
        heights <- c(heights, 1.5*vcfs_per_sample[[i]][j])
      }
      # depth
      if (plot_params$depth) { heights <- c(heights, 7) }
      # breakpoint zoom illustration
      if (split){ heights <- c(heights, 2) }
      # add tracks according to plot params
      if (plot_params$ins) { heights <- c(heights, 4, 4)  }
      if (plot_params$hardclipped) {  heights <- c(heights, 1) }
      if (plot_params$secondary) { heights <- c(heights, 1) }
      if (plot_params$supplementary) { heights <- c(heights, 1) }
      if (plot_params$orphaned) { heights <- c(heights, 1) }
      if (plot_params$inverted) { heights <- c(heights, 1) }
      if (plot_params$samestrand) { heights <- c(heights, 1) }
      # breakpoint zoom axes
      if (split){ heights <- c(heights, 2) }
      # separator
      heights <- c(heights, 1)
    }
    # add in heights for SV_AF tracks
    if (plot_params$svAF) {
      for (i in 1:length(annotations$SV_AF)) {
        heights <- c(heights, 1.5*max(annotations$AF_tracks[[i]]))
      }
    }
    # add in heights for refGene annotation tracks
    if (plot_params$refgene) {
      if (!is.null(annotations$Genes)) {
        heights <- c(heights, 1*nrow(annotations$Genes))
      } else {
        heights <- c(heights, 1)
      }
    }
    # add in bottom x-axis
    heights <- c(heights, 3)
    # add room for legend
    if (plot_params$legend) {
      heights <- c(heights, 6)
    }
    return(heights)
}

# main method
visualise <- function(folder, sample_names, args, outfile, title='') {
  num_samples <- length(sample_names)
  samples <-  lapply(sample_names, function(x) Sample(folder, x))
  vcfs_per_sample <- lapply(samples, function(x) sapply(x$SVtracks, function(y) max(y)))
  xlims <- samples[[1]]$Xlims
  Ins_ylim <- max(sapply(samples, function(x) x$Ins_ylim))
  annotations <- Annotations(folder)
  plot_params <- PlotParams(args)
  heights <- get_plot_layout(plot_params, annotations, num_samples, vcfs_per_sample, samples[[1]]$Split)
  pdf(outfile, title='SVPV Graphics Output', width = 8, height = 0.15* sum(heights), bg = 'white')
  # initialise first layout
  layout(matrix(1:length(heights), length(heights), 1, byrow = TRUE), heights=heights)
  par(mar = c(0.1, 6, 0.1, 2), oma = c(1, 0.1, 0.1, 0.1))
  # plot top x-axis
  add_position_axis(xlims, 3)
  # plot title
  empty_plot(c(0,1))
  par(font = 2)
  text(0.5, 0.5, title, cex = 1.25)
  par(font = 1)
  separator()
  # plot samples
  for (i in 1:num_samples) {
    plot_sample(samples[[i]], plot_params, Ins_ylim)
  }
  # add in heights for SVMAF tracks
  if (plot_params$svAF) {
    for (i in 1:length(annotations$SV_AF)) {
      plot_svs(annotations$SV_AF[[i]], xlims, annotations$AF_tracks[[i]])
    }
  }
  # add in heights for refGene annotation tracks
  if (plot_params$refgene) {  plot_refgenes(annotations$Genes, xlims) }
  # plot bottom x-axis
  add_position_axis(xlims, 1)
  # add legend
  if (plot_params$legend) { add_legend() }
  # add details
  plot_details(samples[[1]]$Bin_size, samples[[1]]$Num_bins)
  graphics.off()
}
# read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
sample_names <- strsplit(as.character(args[1]), ',')[[1]]
folder <- args[2]
outfile <- args[3]
title <- args[4]
visualise(folder, sample_names, args[5:length(args)], outfile, title)
