# parse plot attributes
PlotAttr <- function(folder){
  attr <- read.table(paste0(folder, 'plot_attr.tsv'), header=TRUE, colClasses=c('character', 'numeric', 'numeric'), sep='\t')
  return(list(region=Region(attr$region),
              r_bin_size=attr$r_bin_size,
              r_bin_num=attr$r_bin_num,
              loci=lapply(unlist(strsplit(attr$loci, ',')), function(x) Region(x)),
              l_bin_size=attr$l_bin_size,
              l_bin_num=attr$l_bin_num))
}
# store plot parameters
PlotParams <- function(folder, args) {
  tracks <- list(
    depth = ('-d' %in% args),
    orphaned = ('-or' %in% args),
    inverted = ('-v' %in% args),
    samestrand = ('-ss' %in% args),
    secondary = ('-se' %in% args),
    supplementary = ('-su' %in% args),
    diffmol = ('-dm' %in% args),
    clipped = ('-cl' %in% args),
    ins = ('-i' %in% args),
    refgene = ('-r' %in% args),
    svAF = ('-af' %in% args),
    legend = ('-l' %in% args)
  )
  attr <- PlotAttr(folder)
  if (!is.na(attr$r_bin_size)){
      if (is.na(attr$l_bin_size)) {
        type = 'contiguous' 
      } else { 
        type = 'zoom'
      }
  } else {
    type = 'split'
  }
  if (is.na(attr$l_bin_num)) { num_loci = 0 } else { num_loci = length(attr$loci) }
  return(list(tracks = tracks,
              type = type,
              num_loci = num_loci,
              Attr = attr))
}
# split a region string (eg 'chr1:100-200')
Region <- function(region){
  start <- as.numeric(gsub('^.+:', '', gsub('-.+$', '', region, perl=TRUE), perl=TRUE))
  end <- as.numeric(gsub('^.+-', '', region, perl=TRUE))
  return(list(
    chr=gsub(':.+$', '', region, perl=TRUE),
    start=start,
    end=end,
    xlims = c(start, end)
    ))
}
# parse insert size files
Inserts <- function(folder){
  fwd_ins <- list()
  rvs_ins <- list()
  for (f in list.files(folder)){
    if (grepl('fwd_ins.csv', f)){
      pos <- sub('.fwd_ins.csv', '', basename(f))
      fwd_ins[[pos]] <- scan(paste0(folder, f), sep='\n', what='character', quiet=TRUE)
      fwd_ins[[pos]] <- lapply(fwd_ins[[pos]], function(x) as.numeric(unlist(strsplit(x, ','))))
      rvs_ins[[pos]] <- scan(sub('fwd', 'rvs', paste0(folder, f)), sep='\n', what='character', quiet=TRUE)
      rvs_ins[[pos]] <- lapply(rvs_ins[[pos]], function(x) as.numeric(unlist(strsplit(x, ','))))
    }
  }
  ylim <- estimate_ylim(c(unlist(fwd_ins, use.names=FALSE), unlist(rvs_ins, use.names=FALSE)))
  return(list(fwd=fwd_ins, rvs=rvs_ins, ylim=ylim))
}
# parse depth files
Depths <- function(folder){
  region <- NULL
  loci <- list()
  for (f in list.files(folder)){
    if (grepl('region_depths.tsv', f)){
      region <- read.delim(paste0(folder, f), header = TRUE,  sep = '\t')
    } else if (grepl('.depths.tsv', f)){
      pos <- sub('.depths.tsv', '', basename(f))
      loci[[pos]] <- read.delim(paste0(folder, f), header = TRUE,  sep = '\t')
    }
  }
  return(list(region=region, loci=loci))
}
# parse alignment stats files
AlnStats <- function(folder){
  aln_stats <- list()
  for (f in list.files(folder)){
    if (grepl('.aln_stats.tsv', f)){
      pos <- sub('.aln_stats.tsv', '', basename(f))
      aln_stats[[pos]] <- read.delim(paste0(folder, f), header = TRUE,  sep = '\t')
    }
  }
  return(aln_stats)
}
# parse SV info
VCFs <- function(file){
  if (file.exists(file)) {
    vcfs <- read.delim(file, header = TRUE, sep = '\t')
    vcfs <- split.data.frame(vcfs, vcfs$vcf)
    for (n in names(vcfs)){
      tracks <- get_tracks(vcfs[[n]]$start, vcfs[[n]]$end, vcfs[[n]]$chrom)
      vcfs[[n]] <- list(calls=vcfs[[n]][2:ncol(vcfs[[n]])],
                     tracks = tracks,
                     n_tracks = max(tracks),
                     name = n)
    }
  } else {vcfs = NA}
  return(vcfs)
}
# store data for a given sample
Sample <- function(folder, sample, label = FALSE) {
  my_folder = paste0(folder, sample, '/')
  return(list(
    Name = sample,
    Depths = Depths(my_folder),
    vcfs = VCFs(paste0(my_folder, 'svs.tsv')),
    Ins = Inserts(my_folder),
    aln_stats = AlnStats(my_folder)))
}
# container for annotations
Annotations <- function(folder) {
  #check if refgenes file exists, if so read
  genes_file <- paste0(folder, 'refgene.tsv')
  if (file.exists(genes_file)) {
    genes <- read.delim(genes_file, header=TRUE, sep='\t')
  } else {
    genes = NULL
  }
  return(list(Genes=genes, vcfs=VCFs(paste0(folder, 'SV_AF.tsv'))))
}
# add main plot title
add_title <- function(title, cex=1.25, font=2){
  empty_plot(c(0,1)); par(font = font); text(0.5, 0.5, title, cex =cex); par(font = 1)
}
# add a border around a given plot area
add_border <- function(xlims, ylims, lwd=0.5) {
  rect(xlims[1], ylims[1], xlims[2], ylims[2], lwd=lwd)
}
# horizontal line separator between samples
separator <- function() {
  empty_plot(c(0, 1)); par(xpd=NA); abline(h=0.5, lwd=2, col='gray25'); par(xpd=FALSE)
}
# create an empty plot
empty_plot <- function(xlim,ylim=c(0, 1),type='n',bty='n', xaxt='n', yaxt='n', ylab='', xlab='') {
    plot(1, type=type, ylim=ylim, xlim=xlim, bty=bty, xaxt=xaxt, yaxt=yaxt, ylab=ylab, xlab=xlab)
}
# add axis to a plot
add_position_axis <- function(params, side) {
  if (params$type == 'split'){
    for (i in 1:2){
      units <- get_units(params$Attr$loci[[i]]$end - params$Attr$loci[[i]]$start + 1)
      empty_plot(c(params$Attr$loci[[i]]$start, params$Attr$loci[[i]]$end) / units$val)
      mtext(paste0(params$Attr$loci[[i]]$chr, ' pos (', units$sym, ')'), side=side, line=-1, cex=0.85)
      axis(side=side, line=-3)
    }
  } else {
    units <- get_units(params$Attr$region$end - params$Attr$region$start + 1)
    empty_plot(c(params$Attr$region$start, params$Attr$region$end) / units$val)
    mtext( paste0(params$Attr$region$chr, ' pos (', units$sym, ')'), side=side, line=-1, cex=0.85)
    axis(side=side, line=-3)
  }
}
# plot read depth and mapping quality
plot_depth <- function(params, depths) {
  par(las=1)
  if (params$type != 'split'){
    ylims = c(0, 1.2 * max(depths$region$total - depths$region$mapQ0 - depths$region$mapQltT, na.rm=TRUE))
    empty_plot(params$Attr$region$xlims, ylim=ylims)
    add_border(params$Attr$region$xlims, ylims)
    depth_plotter(depths$region, params$Attr$r_bin_size)
    title(ylab='Depth\n(reads/ bp)', line=2)
    axis(2, tick=TRUE, labels=TRUE, line=0.5)
  } else {
    ylims = c(0, 1.2 * max(sapply(depths$loci, function(x) max(x$total - x$mapQ0 - x$mapQltT, na.rm=TRUE))))
    for (i in 1:2){
      empty_plot(params$Attr$loci[[i]]$xlims, ylim=ylims)
      add_border(params$Attr$loci[[i]]$xlims, ylims)
      depth_plotter(depths$loci[[i]], params$Attr$l_bin_size)
      if (i == 1){ title(ylab='Depth\n(reads/ bp)', line=2); axis(2, tick=TRUE, labels=TRUE, line=0.5) }
    }
  }
}

depth_plotter <- function(depth, bin_size){
  rect(depth$bin, c(0), (depth$bin + bin_size), depth$total, col='#74C476')
  # add depth$mapQ0 to depth
  rect(depth$bin, depth$total, (depth$bin + bin_size), (depth$total - depth$mapQ0), col='white')
  # add depth$mapQltT to depth
  rect(depth$bin, (depth$total - depth$mapQ0), (depth$bin + bin_size), depth$total - (depth$mapQ0 + depth$mapQltT), col='khaki1')
}

# add the legend
add_legend <- function() {
  # create a plot with room for four legends: depth, inserts, mapping stats, svtype/freq
  empty_plot(c(0, 4), ylim=c(0, 1))
  add_border(c(0, 1), c(0, 1))
  add_border(c(1, 2), c(0, 1))
  add_border(c(2, 3), c(0, 1))
  add_border(c(3, 4), c(0, 1))
  par(font=2)
  text(c(0.5, 1.5, 2.5, 3.5), 1,  pos=1, labels=c("Read MapQ", "Inferred Insert Size", "Mapping Stats", "SV Allele Frequency"))
  par(font=1)
  # constants
  bottom = 0.20
  top = bottom + 0.20
  # depth legend
  rect(c((1 / 6 - 0.1), (3 / 6 - 0.1), (5 / 6 - 0.1)), bottom, c((1 / 6 + 0.1), (3 / 6 + 0.1), (5 / 6 + 0.1)), top,  col=c('#74C476', 'khaki1', 'gray95'))
  text( c((1 / 6), (3 / 6), (5 / 6)), c(top + 0.1), pos=3, labels=c(">= 30", "< 30", "= 0") )
  # inferred insert size legend
  text(1.5, top + 0.25, labels=c("Proportion in position x\nwith mapping distance y"))
  rect(1.15 + 0.7 / 10 * (0:9), bottom, 1.15 + 0.7 / 10 * (1:10), top, col=insert_size_pallete(10)[(1:10)])
  text(c(1.18, 1.82), bottom - 0.08, as.character(c(0, 1)))
  # Mapping stats legend 
  text(2.5, top + 0.1, pos=3, labels=c("proportion of reads in position x"))
  rect(2.15 + 0.7 / 10 * (0:9), bottom, 2.15 + 0.7 / 10 * (1:10), top, col=aln_stats_pallete(10)[(1:10)])
  text(c(2.18, 2.82), bottom - 0.08, as.character(c(0, 1)))
  # SV AF legend
  sv_types <- c("DEL", "DUP", "CNV", "INV")
  height = 0.13
  top = bottom + 0.15 * length(sv_types)
  par(family='mono', font=2)
  for (i in 1:length(sv_types)){
    rect(3.20 + 0.7 / 10 * (0:9), bottom + (i-1) * height, 3.20 + 0.7 / 10 * (1:10), bottom + (i) * height, col=sapply(1:10, function(x) get_sv_col(sv_types[i], x/10)))
    text(3.20, bottom + (i-0.5) * height, sv_types[i], pos=2)
  }
  par(family='sans', font=1)
  text(c(3.23, 3.87), bottom - 0.08, as.character(c(0, 1)))
}
aln_stats_pallete <- function(n){
  return(colorRampPalette(c("gray95", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603", "#7F2704"))(n))
}
insert_size_pallete <- function(n){
  return(colorRampPalette(c("gray95", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B"))(n))
}
plot_aln_stats <- function(params, aln_stats){
  if (params$tracks$clipped){ 
    aln_stats_plotter( aln_stats$reads, aln_stats$clipped, 'clipped') 
    }
  if (params$tracks$secondary){ 
    aln_stats_plotter( aln_stats$reads, aln_stats$secondary, 'secondary') 
    }
  if (params$tracks$supplementary){ 
    aln_stats_plotter( aln_stats$reads, aln_stats$supplementary, 'supplementary') 
    }
  if (params$tracks$diffmol){ 
    aln_stats_plotter( aln_stats$reads, aln_stats$supplementary, 'diffmol') 
    }
  if (params$tracks$orphaned){ 
    aln_stats_plotter( aln_stats$reads, aln_stats$orphaned, 'orphaned') 
    }
  if (params$tracks$inverted){ 
    aln_stats_plotter( aln_stats$reads, aln_stats$inverted, 'inverted') 
    }
  if (params$tracks$samestrand){ 
    aln_stats_plotter( aln_stats$reads, aln_stats$samestrand, 'samestrand') 
    }
}
# plot alignment stats
aln_stats_plotter <- function(total, numerator, label, split, spacer=2) {
    if (split) {end=length(total)+4*spacer} else { end=length(total)}
    empty_plot(c(0, end))
    par(las=1)
    mtext(label, side=2, line=-1,  cex=0.75)
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
    return(list(val=1, sym='bp'))
  } else if (num_bp < 1000000) {
    return(list(val=1000, sym='kbp'))
  } else if (num_bp < 1000000000) {
    return(list(val=1000000, sym='Mbp'))
  } else {
    return(list(val=1000000000, sym='Gbp'))
  }
}
# return the 1.1 * lowest insert size in the highest top 15% - heuristic for avoiding outliers
estimate_ylim <- function(ins) {
  sorted <- sort(ins)
  return(1.1 * sorted[floor(length(sorted) * 0.85)])
}
# plots the inserts in specifiend interval
plot_binned_inserts <- function(binned_inserts, num_y_bins, split, spacer=2){
  if (split) {end=nrow(binned_inserts) + 4*spacer} else {end=nrow(binned_inserts)}
  empty_plot(c(0,end), ylab='' ,ylim=c(0, num_y_bins + 1))
  if (split){
    for (i in 1:(num_y_bins + 1)) {
      rect(spacer:(spacer+split-2), i - 1, (spacer+1):(spacer+split-1), i, col=insert_size_pallete(25)[24 * binned_inserts[1:(split-1), i] + 1],  border=NA)
      rect((3*spacer+split-1):(3*spacer+nrow(binned_inserts)-1), i - 1, (3*spacer+split):(3*spacer+nrow(binned_inserts)), i, col=insert_size_pallete(25)[24 * binned_inserts[(split:nrow(binned_inserts)), i] + 1],  border=NA)
      add_border(c(spacer, spacer+split-1), c(0, num_y_bins))
      add_border(c(spacer, spacer+split-1), c(num_y_bins, num_y_bins + 1))
      add_border(c(3*spacer+split-1, 3*spacer+nrow(binned_inserts)), c(0, num_y_bins))
      add_border(c(3*spacer+split-1, 3*spacer+nrow(binned_inserts)), c(num_y_bins, num_y_bins + 1))
    }
  } else {
    for (i in 1:(num_y_bins + 1)) {
      rect((1:nrow(binned_inserts))-1, i - 1, (1:nrow(binned_inserts)), i, col=insert_size_pallete(25)[24 * binned_inserts[(1:nrow(binned_inserts)), i] + 1],  border=NA)
      add_border(c(0, end), c(0, num_y_bins))
      add_border(c(0, end), c(num_y_bins, num_y_bins + 1))
    }
  }
}

plot_insert_sizes <- function(params, ins, ins_ylim){
  # (fwd_ins, rvs_ins, ylim, split, num_y_bins=10) 
  # divide into 10 bins spaced equally between 0 and ylim
  ybin_size <- ylim / num_y_bins
  # create an extra bin to store anythin larger than ylim
  fwd_bins <- matrix(nrow=nrow(fwd_ins), ncol=(num_y_bins + 1))
  rvs_bins <- matrix(nrow=nrow(rvs_ins), ncol=(num_y_bins + 1))
  # get counts for each bin
  for (i in 1:num_y_bins) {
    fwd_bins[, i]=sapply(fwd_ins[,2], function(x) sum((((i - 1) * ybin_size)  <= x) & (x < ((i) * ybin_size))))
    rvs_bins[, i]=sapply(rvs_ins[,2], function(x) sum((((i - 1) * ybin_size)  <= x) & (x < ((i) * ybin_size))))
  }
  fwd_bins[, (num_y_bins + 1)]=sapply(fwd_ins[,2], function(x) sum(ylim <= x))
  rvs_bins[, (num_y_bins + 1)]=sapply(rvs_ins[,2], function(x) sum(ylim <= x))
  # convert to proportions
  fwd_bins=fwd_bins / rowSums(fwd_bins)
  rvs_bins=rvs_bins / rowSums(rvs_bins)
  # organise sensible units for ticks on plot
  units <- get_units(ylim / 2)
  ylim <- ylim / units$val
  mid <- round(ylim / 2)
  interval <- round(mid * 2 / 3)
  ticks_at <- c((mid - interval), mid, (mid + interval))
  # plot binned foward inserts
  par(las=1)
  plot_binned_inserts(fwd_bins, num_y_bins, split)
  axis(2, at=10 * ticks_at / ylim,  labels=as.character(ticks_at),  line=-1)
  title(ylab=paste0('forward\nmapping\ndistance\n (', units$sym, ')'), line=1.5)
  par(xpd=NA)
  text(0, num_y_bins + 0.5, labels =">", cex=0.85, pos=2)
  # plot binned reverse inserts
  plot_binned_inserts(rvs_bins, num_y_bins, split)
  axis( 2, at=10 * ticks_at / ylim,  labels=as.character(ticks_at), line=-1)
  title(ylab=paste0('reverse\nmapping\ndistance\n (', units$sym, ')'), line=1.5)
  text(0, num_y_bins + 0.5, labels =">", cex=0.85, pos=2)
  par(xpd=FALSE)
}
# plot structural variants
#plot_svs(params, sample$svs, AF=FALSE)
sv_plotter <- function(svs, tracks, xlims, AF=TRUE){
  scale=1 / max(tracks)
  par(las=1)
  par(font=2)
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
      rect(start, bottom + spacer, end, top - spacer, col=get_sv_col(svs$svtype[i], as.numeric(svs$AF[i])), border=get_sv_col(svs$svtype[i], 1), lwd=2)
      if (x_prop > 1/5){
        text(0.5 * (start + end), 0.5 * (top + bottom), labels = paste0(svs$svtype[i], ' : AF = ', as.character(round(as.numeric(svs$AF[i]), digits = 3)), ' : ', as.character(round(len/units$val, digits=2)), " ", units$sym))
      } else if (x_prop > 1/10){
        text(0.5 * (start + end), 0.5 * (top + bottom), labels = paste0('AF = ', as.character(round(as.numeric(svs$AF[i]), digits = 3))))
      }
    }
  }
  par(font = 1)
}
# coordinate sv plotting
plot_svs <- function(params, vcf, AF=TRUE) {
  if (params$type != 'split'){
    empty_plot(params$Attr$region$xlims)
    add_border(params$Attr$region$xlims ,c(0,1))
    mtext(vcf$name, side=2, line=0, cex=0.8)
    # text(0.5 * (params$Attr$region$xlims[1] + params$Attr$region$xlims[2]), 0.5, labels="None")
    sv_plotter(vcf$calls, vcf$tracks, params$Attr$region$xlims, AF=AF)
  } else {
    for (i in 1:2){
      empty_plot(params$Attr$loci[[i]]$xlims)
      add_border(params$Attr$loci[[i]]$xlims, c(0,1))
      if (i == 1) mtext(vcf$name, side=2, line=0, cex=0.8)
      idxs = which(vcf$calls$chrom == params$Attr$loci[[i]]$chr)
      if (!length(idxs)) {
        text(0.5 * (params$Attr$loci[[i]]$xlims[1] + params$Attr$loci[[i]]$xlims[2]), 0.5, labels="None")
      } else {
        sv_plotter(vcf$calls[idxs,], vcf$tracks[idxs], params$Attr$loci[[i]]$xlims, AF=AF)
      }
    }
  }
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
plot_sample <- function(sample, params, ins_ylim) {
  add_title(paste("Sample:", sample$Name), cex=1) # plot sample title
  #plot sample SVs
  for (vcf in sample$vcfs) {
    plot_svs(params, vcf, AF=FALSE)
  }
  # plot sample depth
  if (params$tracks$depth) { plot_depth(params, sample$Depths) }
  # plot zoom details
  if (params$type == 'zoom'){ add_zoom_detail(params) }
  # plot sample insert sizes
  if (params$tracks$ins){ plot_insert_sizes(params, sample$Ins, ins_ylim) }
  # plot aln stats
  plot_aln_stats(params, sample$aln_stats)
  # plot zoom details
  if (params$type == 'zoom'){ add_zoom_detail(params,axes=TRUE) }
  # add separator
  separator()
}
# add annotations of details of the plot
plot_details <- function(bin_size, num_bins) {
  mtext( paste("Bin size: ", as.character(bin_size), "    Num bins: ", as.character(num_bins), "    Date: ", as.character(Sys.Date()) ), side = 1, line = 0, adj = 0, cex = 0.65 )
}
# graphical representation of linear transformation between depth plot and zoomed in inserts size and aln stats plots
add_zoom_detail <- function(params, col='black', spacer=2, axes = FALSE){
  #(xlims, start_1, end_1, start_2, end_2, num_bins, col='black', spacer=2, axes = FALSE)
  xlims <- params$Attr$region$xlims
  start_1 <- params$Attr$loci[[1]]$start
  start_2 <- params$Attr$loci[[2]]$start
  end_1 <- params$Attr$loci[[1]]$end
  end_2 <- params$Attr$loci[[2]]$end
  num_bins <- params$Attr$l_bin_num
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
get_tracks <- function(starts, ends, chroms) {
  # assume that starts and ends are of same length and are sorted by lowest start first
  tracks <- vector("integer", length = length(starts))
  tracks[1] = 1
  if (length(starts) >= 2) {
    for (i in 2:length(starts)) {
      for (j in 1:i) {
        overlap = FALSE
        if (j %in% tracks) {
          check = which(tracks %in% j)
          for (k in 1:length(check)) {
            if ((starts[i] <= ends[check[k]]) & (chroms[i] == chroms[check[k]])) {
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
# add refgene tracks to the plot
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
# returns the appropriate matrix and heights for the plot layout
get_plot_layout <- function(params, annotations, num_samples, vcfs_per_sample, max=200) {
  if (params$type != 'split') n_col = 1 else n_col = 2
  # generate vectors heights h, plot number p_n for layout
  # top x-axis
  h <- c(4); p_n <- 1:n_col; idx <- n_col
  # title
  h <- c(h,1); p_n <- c(p_n, rep(p_n[idx]+1, times=n_col)); idx=idx+n_col
  # separator
  h <- c(h, 1); p_n <- c(p_n, rep(p_n[idx]+1, times=n_col)); idx=idx+n_col
  for (i in 1:num_samples){
    # sample title
    h <- c(h, 1.5); p_n <- c(p_n, rep(p_n[idx]+1, times=n_col)); idx=idx+n_col
    # add in SV plots
    for (j in 1:length(vcfs_per_sample[[i]])){
      h <- c(h, 1.5*vcfs_per_sample[[i]][j]); p_n <- c(p_n, (p_n[idx]+1):(p_n[idx]+n_col)); idx=idx+n_col
    }
    # depth
    if (params$tracks$depth) { h <- c(h, 7); p_n <- c(p_n, (p_n[idx]+1):(p_n[idx]+n_col)); idx=idx+n_col }
    # breakpoint zoom illustration
    if (params$type == 'zoom'){ h <- c(h, 2); p_n <- c(p_n, (p_n[idx]+1):(p_n[idx]+n_col)); idx=idx+n_col }
    # add tracks according to plot params
    if (params$tracks$ins) { h <- c(h, 4, 4); p_n <- c(p_n, (p_n[idx]+1):(p_n[idx]+2*n_col)); idx=idx+2*n_col  }
    if (params$tracks$clipped) {  h <- c(h, 1); p_n <- c(p_n, (p_n[idx]+1):(p_n[idx]+n_col)); idx=idx+n_col }
    if (params$tracks$secondary) { h <- c(h, 1); p_n <- c(p_n, (p_n[idx]+1):(p_n[idx]+n_col)); idx=idx+n_col }
    if (params$tracks$supplementary) { h <- c(h, 1); p_n <- c(p_n, (p_n[idx]+1):(p_n[idx]+n_col)); idx=idx+n_col }
    if (params$tracks$diffmol) { h <- c(h, 1); p_n <- c(p_n, (p_n[idx]+1):(p_n[idx]+n_col)); idx=idx+n_col }
    if (params$tracks$orphaned) { h <- c(h, 1); p_n <- c(p_n, (p_n[idx]+1):(p_n[idx]+n_col)); idx=idx+n_col }
    if (params$tracks$inverted) { h <- c(h, 1); p_n <- c(p_n, (p_n[idx]+1):(p_n[idx]+n_col)); idx=idx+n_col }
    if (params$tracks$samestrand) { h <- c(h, 1); p_n <- c(p_n, (p_n[idx]+1):(p_n[idx]+n_col)); idx=idx+n_col }
    # breakpoint zoom axes
    if (params$type == 'zoom'){ h <- c(h, 2); p_n <- c(p_n, (p_n[idx]+1):(p_n[idx]+n_col)); idx=idx+n_col }
    # separator
    h <- c(h, 1); p_n <- c(p_n, rep(p_n[idx]+1, times=n_col)); idx=idx+n_col
  }
  # add in h for SV_AF tracks
  if (params$tracks$svAF) {
    for (vcf in annotations$vcfs) {
      h <- c(h, 1.5*vcf$n_tracks); p_n <- c(p_n, (p_n[idx]+1):(p_n[idx]+n_col)); idx=idx+n_col
    }
  }
  # add in h for refGene annotation tracks
  if (params$tracks$refgene) {
    if (length(annotations$Genes$chrom) > 1) {
      h <- c(h, 1*nrow(annotations$Genes)); p_n <- c(p_n, (p_n[idx]+1):(p_n[idx]+n_col)); idx=idx+n_col
    } else {
      h <- c(h, 1); p_n <- c(p_n, (p_n[idx]+1):(p_n[idx]+n_col)); idx=idx+n_col
    }
  }
  # add in bottom x-axis
  h <- c(h, 4); p_n <- c(p_n, (p_n[idx]+1):(p_n[idx]+n_col)); idx=idx+n_col
  # add room for legend
  if (params$tracks$legend) {
    h <- c(h, 6); p_n <- c(p_n, rep(p_n[idx]+1, times=n_col)); idx=idx+n_col
  }
  mat <- cbind((p_n[idx]+1):(p_n[idx]+length(h)), matrix(p_n, length(h), n_col, byrow = TRUE))
  widths <- c(1, rep(8/n_col, times=n_col))
  return(list(mat=mat, heights=h, widths=widths))
}
# main method
visualise <- function(folder, sample_names, plot_args, outfile, title='') {
  params <- PlotParams(folder, plot_args)
  num_samples <- length(sample_names)
  samples <-  lapply(sample_names, function(x) Sample(folder, x))
  vcfs_per_sample <- lapply(samples, function(x) sapply(x$vcfs, function(y) y$n_tracks))
  ins_ylim <- max(sapply(samples, function(x) x$Ins$ylim))
  annotations <- Annotations(folder)
  lay_out <- get_plot_layout(params, annotations, num_samples, vcfs_per_sample)
  # pdf(outfile, title='SVPV Graphics Output', width = 8, height = 0.15* sum(lay$heights), bg = 'white')
  layout(lay_out$mat, heights=lay_out$heights, widths=lay_out$widths)
  par(mai = c(0.01, 0.01, 0.01, 0.01), omi = c(0,0,0,0))
  add_position_axis(params, 3) # top x axis
  add_title(title)
  separator()
  # plot samples
  for (i in 1:num_samples) {
    plot_sample(samples[[i]], params, ins_ylim)
  }
  # add in heights for SVAF tracks
  if (params$tracks$svAF) {
    for (i in 1:length(annotations$SV_AF)) {
      plot_svs(annotations$SV_AF[[i]], xlims, annotations$AF_tracks[[i]])
    }
  }
  # add in heights for refGene annotation tracks
  if (params$tracks$refgene) {  plot_refgenes(annotations$Genes, xlims) }
  # plot bottom x-axis
  add_position_axis(xlims, 1)
  # add legend
  if (params$tracks$legend) { add_legend() }
  # add details
  plot_details(samples[[1]]$Bin_size, samples[[1]]$Num_bins)
  graphics.off()
}
# read command-line arguments
#args <- commandArgs(trailingOnly = TRUE)
args <- c('NA12877_S1,NA12878_S1,NA12884_S1', '/home/jacmun/R/SVPV_debug/INV/chr1_247127925/', '/home/jacmun/R/SVPV_debug/INV/chr1_247127925/cdfbdc4880.pdf', "DEL at chr1:4064686-4066276", '-d', '-or', '-v', '-ss', '-cl', '-i', '-r', '-af', '-l', '-dm')
sample_names <- strsplit(as.character(args[1]), ',')[[1]]
folder <- args[2]
outfile <- args[3]
title <- args[4]
plot_args <- args[5:length(args)]
visualise(folder, sample_names, plot_args, outfile, title)
