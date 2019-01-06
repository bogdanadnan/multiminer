#!/usr/bin/Rscript

require(ggplot2)
require(scales)

NS_PER_SEC <- 1000 * 1000 * 1000

theme_set(theme_gray(base_size = 10))

GRAPH_WIDTH_CM <- 10
GRAPH_HEIGHT_CM <- 3

save_graph <- function(name, plot) {
  ggsave(paste0(name, '.pdf'), plot, width = GRAPH_WIDTH_CM, height = GRAPH_HEIGHT_CM, scale = 2.3, units = 'cm')
}

args <- commandArgs(trailingOnly = TRUE)

bench_id <- args[1]
commits <- args[2:length(args)]

### Load data:
data <- data.frame()

for (bench_mode in c('t_cost', 'm_cost')) {
  for (commit in commits) {
    file <- paste0('bench-', bench_mode, '-', bench_id, '-', commit, '.csv')
    file_data <- read.csv(file)
    data <- rbind(data, data.frame(Bench.mode=bench_mode, Commit=commit, Mode=file_data$mode, Kernel.mode=file_data$kernel, Version=file_data$version, Type=file_data$type, Precompute=file_data$precompute, t_cost=file_data$t_cost, m_cost=file_data$m_cost, lanes=file_data$lanes, ns_per_hash=file_data$ns_per_hash))
  }
}

### Compute additional columns:
data$hashes_per_second <- NS_PER_SEC / data$ns_per_hash
data$blocks_per_second <- data$m_cost * data$t_cost * data$hashes_per_second

make_plots_commits <- function(mode, kernel, type, precompute) {
  data_b <- data[data$Mode == mode & data$Kernel.mode == kernel & data$Version == 'v1.3' & data$Type == paste0('Argon2', type) & data$Precompute == precompute,]
  
  prefix <- paste0('plot-commits-', bench_id, '-', mode, '-', kernel, '-argon2', type)
  if (precompute == 'yes') {
      prefix <- paste0(prefix, '-precompute')
  }
  
  data_t_cost <- data_b[data_b$Bench.mode == 't_cost',]
  data_m_cost <- data_b[data_b$Bench.mode == 'm_cost',]
  
  if (length(data_t_cost$hashes_per_second) != c(0)) {
    save_graph(paste0(prefix, '-t_cost'),
               ggplot(data_t_cost, aes(x=t_cost, y=blocks_per_second, group=Commit, colour=Commit)) +
                 geom_line() +
                 scale_y_continuous(labels=comma, limits=c(0, NA)) +
                 facet_grid(~lanes, labeller=label_both) +
                 xlab('t_cost') + ylab('Blocks per second'))
  }

  if (length(data_m_cost$hashes_per_second) != c(0)) {
    save_graph(paste0(prefix, '-m_cost'),
               ggplot(data_m_cost, aes(x=m_cost, y=blocks_per_second, group=Commit, colour=Commit)) +
                 geom_line() +
                 scale_x_continuous(trans="log2", labels=trans_format("log2", math_format(2^.x))) +
                 scale_y_continuous(labels=comma, limits=c(0, NA)) +
                 facet_grid(~lanes, labeller=label_both) +
                 xlab('m_cost (log scale)') + ylab('Blocks per second'))
  }
}

make_plots_types <- function(commit, mode, kernel) {
  data_b <- data[data$Commit == commit & data$Mode == mode & data$Kernel.mode == kernel & data$Version == 'v1.3',]
  if (length(data_b$hashes_per_second) != c(0)) {
    data_b$Variant <- paste0(data_b$Type, '-', data_b$Precompute)
  }
  
  prefix <- paste0('plot-types-', bench_id, '-', commit, '-', mode, '-', kernel)
  
  data_t_cost <- data_b[data_b$Bench.mode == 't_cost',]
  data_m_cost <- data_b[data_b$Bench.mode == 'm_cost',]
  
  if (length(data_t_cost$hashes_per_second) != c(0)) {
    save_graph(paste0(prefix, '-t_cost'),
               ggplot(data_t_cost, aes(x=t_cost, y=blocks_per_second, group=Variant, colour=Variant)) +
                 geom_line() +
                 scale_y_continuous(labels=comma, limits=c(0, NA)) +
                 facet_grid(~lanes, labeller=label_both) +
                 xlab('t_cost') + ylab('Blocks per second'))
  }

  if (length(data_m_cost$hashes_per_second) != c(0)) {
    save_graph(paste0(prefix, '-m_cost'),
               ggplot(data_m_cost, aes(x=m_cost, y=blocks_per_second, group=Variant, colour=Variant)) +
                 geom_line() +
                 scale_x_continuous(trans="log2", labels=trans_format("log2", math_format(2^.x))) +
                 scale_y_continuous(labels=comma, limits=c(0, NA)) +
                 facet_grid(~lanes, labeller=label_both) +
                 xlab('m_cost (log scale)') + ylab('Blocks per second'))
  }
}

make_plots_versions <-  function(commit, mode, kernel, type, precompute) {
  data_b <- data[data$Commit == commit & data$Mode == mode & data$Kernel.mode == kernel & data$Type == paste0('Argon2', type) & data$Precompute == precompute,]
  
  prefix <- paste0('plot-versions-', bench_id, '-', commit, '-', mode, '-', kernel, '-argon2', type)
  if (precompute == 'yes') {
    prefix <- paste0(prefix, '-precompute')
  }
  
  data_t_cost <- data_b[data_b$Bench.mode == 't_cost',]
  data_m_cost <- data_b[data_b$Bench.mode == 'm_cost',]
  
  if (length(data_t_cost$hashes_per_second) != c(0)) {
    save_graph(paste0(prefix, '-t_cost'),
               ggplot(data_t_cost, aes(x=t_cost, y=blocks_per_second, group=Version, colour=Version)) +
                 geom_line() +
                 scale_y_continuous(labels=comma, limits=c(0, NA)) +
                 facet_grid(~lanes, labeller=label_both) +
                 xlab('t_cost') + ylab('Blocks per second'))
  }

  if (length(data_m_cost$hashes_per_second) != c(0)) {
    save_graph(paste0(prefix, '-m_cost'),
               ggplot(data_m_cost, aes(x=m_cost, y=blocks_per_second, group=Version, colour=Version)) +
                 geom_line() +
                 scale_x_continuous(trans="log2", labels=trans_format("log2", math_format(2^.x))) +
                 scale_y_continuous(labels=comma, limits=c(0, NA)) +
                 facet_grid(~lanes, labeller=label_both) +
                 xlab('m_cost (log scale)') + ylab('Blocks per second'))
  }
}

make_plots_kernels <-  function(commit, mode, version, type, precompute) {
  data_b <- data[data$Commit == commit & data$Mode == mode & data$Version == paste0('v', version) & data$Type == paste0('Argon2', type) & data$Precompute == precompute,]
  
  prefix <- paste0('plot-kernels-', bench_id, '-', commit, '-', mode, '-v', version, '-argon2', type)
  if (precompute == 'yes') {
    prefix <- paste0(prefix, '-precompute')
  }
  
  data_t_cost <- data_b[data_b$Bench.mode == 't_cost',]
  data_m_cost <- data_b[data_b$Bench.mode == 'm_cost',]
  
  if (length(data_t_cost$hashes_per_second) != c(0)) {
    save_graph(paste0(prefix, '-t_cost'),
               ggplot(data_t_cost, aes(x=t_cost, y=blocks_per_second, group=Kernel.mode, colour=Kernel.mode)) +
                 geom_line() +
                 scale_y_continuous(labels=comma, limits=c(0, NA)) +
                 facet_grid(~lanes, labeller=label_both) +
                 xlab('t_cost') + ylab('Blocks per second'))
  }

  if (length(data_m_cost$hashes_per_second) != c(0)) {
    save_graph(paste0(prefix, '-m_cost'),
               ggplot(data_m_cost, aes(x=m_cost, y=blocks_per_second, group=Kernel.mode, colour=Kernel.mode)) +
                 geom_line() +
                 scale_x_continuous(trans="log2", labels=trans_format("log2", math_format(2^.x))) +
                 scale_y_continuous(labels=comma, limits=c(0, NA)) +
                 facet_grid(~lanes, labeller=label_both) +
                 xlab('m_cost (log scale)') + ylab('Blocks per second'))
  }
}

make_plots_modes <-  function(commit, kernel, version, type, precompute) {
  data_b <- data[data$Commit == commit & data$Kernel.mode == kernel & data$Version == paste0('v', version) & data$Type == paste0('Argon2', type) & data$Precompute == precompute,]
  
  prefix <- paste0('plot-modes-', bench_id, '-', commit, '-', kernel, '-v', version, '-argon2', type)
  if (precompute == 'yes') {
    prefix <- paste0(prefix, '-precompute')
  }
  
  data_t_cost <- data_b[data_b$Bench.mode == 't_cost',]
  data_m_cost <- data_b[data_b$Bench.mode == 'm_cost',]
  
  if (length(data_t_cost$hashes_per_second) != c(0)) {
    save_graph(paste0(prefix, '-t_cost'),
               ggplot(data_t_cost, aes(x=t_cost, y=blocks_per_second, group=Mode, colour=Mode)) +
                 geom_line() +
                 scale_y_continuous(labels=comma, limits=c(0, NA)) +
                 facet_grid(~lanes, labeller=label_both) +
                 xlab('t_cost') + ylab('Blocks per second'))
  }

  if (length(data_m_cost$hashes_per_second) != c(0)) {
    save_graph(paste0(prefix, '-m_cost'),
               ggplot(data_m_cost, aes(x=m_cost, y=blocks_per_second, group=Mode, colour=Mode)) +
                 geom_line() +
                 scale_x_continuous(trans="log2", labels=trans_format("log2", math_format(2^.x))) +
                 scale_y_continuous(labels=comma, limits=c(0, NA)) +
                 facet_grid(~lanes, labeller=label_both) +
                 xlab('m_cost (log scale)') + ylab('Blocks per second'))
  }
}

# Compare commits:
for (mode in c('opencl', 'cuda', 'cpu')) {
  for (kernel in c('by-segment', 'oneshot')) {
    for (type in c('i', 'd', 'id')) {
      if (type == 'd') {
        precomputes <- c('no')
      } else {
        precomputes <- c('no', 'yes')
      }
      for (precompute in precomputes) {
        make_plots_commits(mode, kernel, type, precompute)
      }
    }
  }
}

# Compare types:
for (commit in commits) {
  for (mode in c('opencl', 'cuda', 'cpu')) {
    for (kernel in c('by-segment', 'oneshot')) {
      make_plots_types(commit, mode, kernel)
    }
  }
}

# Compare versions:
for (commit in commits) {
  for (mode in c('opencl', 'cuda', 'cpu')) {
    for (kernel in c('by-segment', 'oneshot')) {
      for (type in c('i', 'd', 'id')) {
        if (type == 'd') {
          precomputes <- c('no')
        } else {
          precomputes <- c('no', 'yes')
        }
        for (precompute in precomputes) {
          make_plots_versions(commit, mode, kernel, type, precompute)
        }
      }
    }
  }
}

# Compare kernel types:
for (commit in commits) {
  for (mode in c('opencl', 'cuda', 'cpu')) {
    for (version in c('1.0', '1.3')) {
      for (type in c('i', 'd', 'id')) {
        if (type == 'd') {
          precomputes <- c('no')
        } else {
          precomputes <- c('no', 'yes')
        }
        for (precompute in precomputes) {
          make_plots_kernels(commit, mode, version, type, precompute)
        }
      }
    }
  }
}

# Compare modes:
for (commit in commits) {
  for (kernel in c('by-segment', 'oneshot')) {
    for (version in c('1.0', '1.3')) {
      for (type in c('i', 'd', 'id')) {
        if (type == 'd') {
          precomputes <- c('no')
        } else {
          precomputes <- c('no', 'yes')
        }
        for (precompute in precomputes) {
          make_plots_modes(commit, kernel, version, type, precompute)
        }
      }
    }
  }
}
