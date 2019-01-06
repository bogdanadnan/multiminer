#!/usr/bin/Rscript

require(ggplot2)
require(scales)

NS_PER_MS <- 1000 * 1000
NS_PER_SEC <- 1000 * 1000 * 1000

theme_set(theme_gray(base_size = 10))

GRAPH_WIDTH_CM <- 12
GRAPH_HEIGHT_CM <- 6

save_graph <- function(name, plot) {
  ggsave(paste0(name, '.pdf'), plot, width = GRAPH_WIDTH_CM, height = GRAPH_HEIGHT_CM, scale = 2.3, units = 'cm')
}

machines <- read.csv('machines.csv', row.names=NULL)

machines$id <- paste0(machines$type, '-', machines$id)

data <- data.frame()

for (i in 1:nrow(machines)) {
  name <- paste(machines[i,]$name)
  type <- paste(machines[i,]$type)
  id <- paste(machines[i,]$id)
  power <- paste(machines[i,]$power)
  if (type == 'cpu') {
    branch <- 'master'
    mode <- 'cpu'
  } else {
    branch <- 'warp-shuffle'
    mode <- 'cuda'
  }
  
  file <- paste0('bench-', id, '-', branch, '.csv')
  if (file.exists(file)) {
    data_file <- read.csv(file, row.names=NULL)
    data_file <- data_file[data_file$mode == mode & data_file$kernel == 'by-segment',]
    if (type == 'gpu') {
      data_file <- data_file[data_file$precompute == 'yes' | data_file$type == 'Argon2d',]
    }
    
    data <- rbind(data, data.frame(Machine=id, Name=name, HW=toupper(type), Version=data_file$version, Type=data_file$type, t_cost=data_file$t_cost, m_cost=data_file$m_cost, lanes=data_file$lanes, ns_per_hash=data_file$ns_per_hash, power=power))
  }
}

data$hashes_per_second <- NS_PER_SEC / data$ns_per_hash
data$blocks_per_second <- data$m_cost * data$t_cost * data$hashes_per_second

to_numbers <- function(x) {
  as.numeric(levels(x))[x]
}

data$hashes_per_joule <- data$hashes_per_second / to_numbers(data$power)
data$blocks_per_joule <- data$m_cost * data$t_cost * data$hashes_per_joule

make_plots_bps <- function(type) {
  data_b <- data[data$Version == 'v1.3' & data$Type == paste0('Argon2', type),]
  
  if (length(data_b$hashes_per_second) != c(0)) {
    data_b_f_t_cost <- data_b$t_cost %in% c(1, 2, 4, 8, 16)
    data_b_f_m_cost <- data_b$m_cost %in% c(1024, 4096, 16384, 65536, 262144)
    
    prefix <- paste0('plot-bps-argon2', type)

    save_graph(paste0(prefix, '-t_cost'),
               ggplot(data_b[data_b_f_m_cost,], aes(x=t_cost, y=blocks_per_second, colour=Name)) +
                 geom_line() +
                 scale_y_continuous(labels=comma) +
                 facet_grid(lanes~m_cost, labeller=label_both) +
                 xlab('t_cost') + ylab('Blocks per second'))
    
    save_graph(paste0(prefix, '-m_cost'),
               ggplot(data_b[data_b_f_t_cost,], aes(x=m_cost, y=blocks_per_second, colour=Name)) +
                 geom_line() +
                 scale_x_continuous(trans="log2", labels=trans_format("log2", math_format(2^.x))) +
                 scale_y_continuous(labels=comma) +
                 facet_grid(t_cost~lanes, labeller=label_both) +
                 xlab('m_cost (log scale)') + ylab('Blocks per second'))
    
    save_graph(paste0(prefix, '-lanes'),
               ggplot(data_b[data_b_f_t_cost & data_b_f_m_cost,], aes(x=lanes, y=blocks_per_second, colour=Name)) +
                 geom_line() +
                 scale_y_continuous(labels=comma) +
                 facet_grid(t_cost~m_cost, labeller=label_both) +
                 xlab('lanes') + ylab('Blocks per second'))
  }
}

make_plots_bpj <- function(type) {
  data_b <- data[data$Version == 'v1.3' & data$Type == paste0('Argon2', type),]
  
  if (length(data_b$hashes_per_second) != c(0)) {
    data_b_f_t_cost <- data_b$t_cost %in% c(1, 2, 4, 8, 16)
    data_b_f_m_cost <- data_b$m_cost %in% c(1024, 4096, 16384, 65536, 262144)
    
    prefix <- paste0('plot-bpj-argon2', type)
    
    save_graph(paste0(prefix, '-t_cost'),
               ggplot(data_b[data_b_f_m_cost,], aes(x=t_cost, y=blocks_per_joule, colour=Name)) +
                 geom_line() +
                 scale_y_continuous(labels=comma) +
                 facet_grid(lanes~m_cost, labeller=label_both) +
                 xlab('t_cost') + ylab('Blocks per joule'))
    
    save_graph(paste0(prefix, '-m_cost'),
               ggplot(data_b[data_b_f_t_cost,], aes(x=m_cost, y=blocks_per_joule, colour=Name)) +
                 geom_line() +
                 scale_x_continuous(trans="log2", labels=trans_format("log2", math_format(2^.x))) +
                 scale_y_continuous(labels=comma) +
                 facet_grid(t_cost~lanes, labeller=label_both) +
                 xlab('m_cost (log scale)') + ylab('Blocks per joule'))
    
    save_graph(paste0(prefix, '-lanes'),
               ggplot(data_b[data_b_f_t_cost & data_b_f_m_cost,], aes(x=lanes, y=blocks_per_joule, colour=Name)) +
                 geom_line() +
                 scale_y_continuous(labels=comma) +
                 facet_grid(t_cost~m_cost, labeller=label_both) +
                 xlab('lanes') + ylab('Blocks per joule'))
  }
}

for (type in c('i', 'd', 'id')) {
  make_plots_bps(type)
  make_plots_bpj(type)
}

