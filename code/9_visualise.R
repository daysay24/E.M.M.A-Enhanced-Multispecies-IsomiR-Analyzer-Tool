library(ggplot2)
library(patchwork)
library(dplyr)
library(ggrepel)
library(tidyr)

# Passed arguments
args <- commandArgs(trailingOnly = TRUE)

# Create graph folder 
dir.create(file.path(args[[1]]))

# Visualise graph 1
graph.data.1 <- read.csv(paste(args[[2]], 'graph_1_data.csv', sep=''))
names(graph.data.1)[names(graph.data.1) == "type"] <- "Type"

# Customize x axis labels to show ratio 
graph.data.ratio.1<- graph.data.1 %>% 
  select(Type, rpm, group) %>%
  pivot_wider(names_from=Type, values_from = rpm) %>% 
  mutate(ratio=IsomiR/Canonical) 
x.labels.1 <- c()
for (row in 1:nrow(graph.data.ratio.1)) {
  r = graph.data.ratio.1[row, ]
  x.labels.1 <- append(x.labels.1, paste(r['group'], paste("1", round(r['ratio'], digits=2), sep=":"), sep="\n"))
}
graph.data.1$group <- factor(graph.data.1$group, levels = unique(graph.data.ratio.1$group))

graph.1 <- ggplot(data = graph.data.1, aes(x = group, y = rpm, group = Type, color = Type)) +
  geom_line() + 
  geom_point() + 
  scale_x_discrete(labels=x.labels.1) + 
  xlab("Group\nCanonical:IsomiR ratio") +
  theme_bw()+
  theme(legend.title = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 9), 
        axis.text = element_text(size = 9), 
        axis.title = element_text(size = 10), 
        plot.margin = unit(c(8.5,8.5,8.5,8.5), "pt"))
graph.1 <- graph.1 + ggplot(data = graph.data.1, aes(x = group, y = relative_abundance, group = Type, color = Type)) +
  geom_line() + 
  geom_point() + 
  ylab("%") +
  xlab("Group") +
  theme_bw() +
  theme(legend.title = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 9), 
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        plot.margin = unit(c(8.5,8.5,8.5,8.5), "pt"))
graph.1 <- graph.1 + 
  plot_layout(ncol = 2, guides = "collect") 
ggsave(file=paste(args[[1]], 'graph_1.pdf', sep=''), plot = graph.1, height = 5, width = 10)

# Visualise graph 2
graph.data.2 <- read.csv(paste(args[[2]], 'graph_2_data.csv', sep=''))
names(graph.data.2)[names(graph.data.2) == "grouped_type"] <- "Type"
grouped.graph.data.2 <- graph.data.2 %>% group_by(group) %>% group_data() 
grouped.graph.data.2 <- as.data.frame(grouped.graph.data.2)
graph.2 <- NULL
group.scaling <- length(unique(graph.data.2[, 'group'])) / 2 

for (row in 1:nrow(grouped.graph.data.2)) {
  group.name = grouped.graph.data.2[row,'group']
  indexes = unlist(grouped.graph.data.2[row,'.rows'])
  sub.graph.data.2 <- graph.data.2[indexes,]
  
  if (is.null(graph.2)) {
    graph.2 <- ggplot(sub.graph.data.2, aes(x="", y=rpm, fill=Type)) +
      geom_bar(stat="identity", width=1, color="white") +
      coord_polar("y", start=0) +
      theme_void() +
      geom_label_repel(aes(label = paste(round(unique_tag/sum(unique_tag)*100,1), "%")),
                       color="white", direction = "y", position = position_stack(vjust = .5), size = 3.5, show.legend = FALSE)+
      ggtitle(paste(group.name, "isomiRNAs (rpm)", sep=" "))+
      theme(plot.title = element_text(hjust = 0.5, size = 12.5, face = "bold"), 
            legend.title = element_text(size = 10, face = "bold"), 
            legend.text = element_text(size = 9))
  } else {
    graph.2 <- graph.2 + ggplot(sub.graph.data.2, aes(x="", y=rpm, fill=Type)) +
      geom_bar(stat="identity", width=1, color="white") +
      coord_polar("y", start=0) +
      theme_void() +
      geom_label_repel(aes(label = paste(round(unique_tag/sum(unique_tag)*100,1), "%")),
                       color="white", position = position_stack(vjust = .5), size = 3.5, show.legend = FALSE)+
      ggtitle(paste(group.name, "isomiRNAs (rpm)", sep=" "))+
      theme(plot.title = element_text(hjust = 0.5, size = 12.5, face = "bold"), 
            legend.title = element_text(size = 10, face = "bold"),
            legend.text = element_text(size = 9))
  }
  
  graph.2 <- graph.2 + ggplot(sub.graph.data.2, aes(x="", y=unique_tag, fill=Type)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() +
    geom_label_repel(aes(label = paste(round(unique_tag/sum(unique_tag)*100,1), "%")),
               color="white", position = position_stack(vjust = .5), size = 3.5, show.legend = FALSE)+
    ggtitle(paste(group.name, "isomiRNAs (unique tag)", sep=" "))+
    theme(plot.title = element_text(hjust = 0.5, size = 12.5, face = "bold"),
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 9))
}
graph.2 <- graph.2 + 
  plot_layout(ncol = 2, guides = "collect")
ggsave(file=paste(args[[1]], 'graph_2.pdf', sep=''), plot = graph.2, width = 10, height = 10 * group.scaling)

# Visualise graph 3
graph.data.3 <- read.csv(paste(args[[2]], 'graph_3_data.csv', sep=''))
names(graph.data.3)[names(graph.data.3) == "type_nt"] <- "Type"
ends <- list("All isomiR", "3'isomiR", "5'isomiR")
graph.3 <- NULL
group.scaling <- length(unique(graph.data.3[, 'group'])) / 2 
type_nt.scaling <- length(unique(graph.data.3[, 'Type'])) / 20 
for (end in ends) {
  sub.graph.data.3 <- NULL
  if (end == 'All isomiR') {
    sub.graph.data.3 <- graph.data.3
  } else {
    sub.graph.data.3 <- graph.data.3[which(graph.data.3$grouped_type == end), ] 
  }
  
  if (is.null(graph.3)) {
    graph.3 <- ggplot(data=sub.graph.data.3, aes(x=group, y=rpm, fill=Type)) + 
      ggtitle(paste(end, "(rpm)", sep=" ")) + 
      geom_bar(position = "fill", stat="identity") + 
      ylab("") +  
      xlab("") +
      theme_bw()+
      theme(axis.ticks.x=element_blank(), 
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_blank(),
            plot.title = element_text(size = 25, face = "bold"), 
            axis.text = element_text(size = 18),
            legend.title = element_text(size = 20, face = "bold"),
            legend.text = element_text(size = 18),
            plot.margin = unit(c(17,17,17,17), "pt")) + 
      scale_y_continuous(labels = scales::percent) 
  } else {
    graph.3 <- graph.3 + ggplot(data=sub.graph.data.3, aes(x=group, y=rpm, fill=Type)) + 
      ggtitle(paste(end, "(rpm)", sep=" ")) +
      geom_bar(position = "fill", stat="identity") + 
      ylab("") +  
      xlab("") +
      theme_bw()+
      theme(axis.ticks.x=element_blank(), 
            panel.grid.major.x = element_blank(), 
            panel.grid.major.y = element_blank(),
            plot.title = element_text(size = 25, face = "bold"),
            axis.text = element_text(size = 18),
            legend.title = element_text(size = 20, face = "bold"),
            legend.text = element_text(size = 18),
            plot.margin = unit(c(17,17,17,17), "pt")) + 
      scale_y_continuous(labels = scales::percent) 
  }
  graph.3 <- graph.3 + ggplot(data=sub.graph.data.3, aes(x=group, y=unique_tag, fill=Type)) + 
    ggtitle(paste(end, "(unique tag)", sep=" ")) +
    geom_bar(position = "fill", stat="identity") + 
    ylab("") +  
    xlab("") +
    theme_bw() +
    theme(axis.ticks.x=element_blank(), 
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(size = 25, face = "bold"),
          axis.text = element_text(size = 18),
          legend.title = element_text(size = 20, face = "bold"),
          legend.text = element_text(size = 18),
          plot.margin = unit(c(17,17,17,17), "pt")) + 
    scale_y_continuous(labels = scales::percent) 
}
graph.3 <- graph.3 + 
  plot_layout(ncol = 2, guides = "auto") 
ggsave(file=paste(args[[1]], 'graph_3.pdf', sep=''), plot = graph.3, width = 20 * group.scaling, height = 20 * type_nt.scaling, limitsize = FALSE)

# Visualise graph 4
graph.data.4 <- read.csv(paste(args[[2]], 'graph_4_data.csv', sep=''))
names(graph.data.4)[names(graph.data.4) == "templated"] <- "Templated"
names(graph.data.4)[names(graph.data.4) == "group"] <- "Group"
names(graph.data.4)[names(graph.data.4) == "position"] <- "Position"
grouped.graph.data.4 <- graph.data.4 %>% group_by(Group) %>% group_data() 
grouped.graph.data.4 <- as.data.frame(grouped.graph.data.4)
graph.4 <- NULL

group.scaling <- length(unique(graph.data.4[, 'Group'])) / 2 
type_ext.scaling <- length(unique(graph.data.4[, 'Position'])) / 10 

for (row in 1:nrow(grouped.graph.data.4)) {
  group.name = grouped.graph.data.4[row,'Group']
  indexes = unlist(grouped.graph.data.4[row,'.rows'])
  sub.graph.data.4 <- graph.data.4[indexes,]
  
  sub.graph.data.4$Position=factor(sub.graph.data.4$Position, levels = unique(sub.graph.data.4$Position))
  
  if (is.null(graph.4)) {
    graph.4 <- ggplot(sub.graph.data.4, aes(fill=Templated, y=count, x=Position)) + 
      ggtitle(paste(group.name, "(%)", sep=" ")) +
      geom_bar(position="fill", stat="identity") +
      xlab("Position") +
      ylab("%")+
      theme_bw() + 
      theme(plot.title = element_text(size = 12.5, face = "bold"), 
            legend.title = element_text(size = 10, face = "bold"), 
            legend.text = element_text(size = 9), 
            axis.text = element_text(size = 9), 
            axis.title = element_text(size = 10), 
            plot.margin = unit(c(8.5,8.5,8.5,8.5), "pt")) + 
      scale_y_continuous(labels = scales::percent) 
  } else {
    graph.4 <- graph.4 + ggplot(sub.graph.data.4, aes(fill=Templated, y=count, x=Position)) + 
      ggtitle(paste(group.name, "(%)", sep=" ")) +
      geom_bar(position="fill", stat="identity") +
      xlab("Position") +
      ylab("%") +
      theme_bw() + 
      theme(plot.title = element_text(size = 12.5, face = "bold"), 
            legend.title = element_text(size = 10, face = "bold"), 
            legend.text = element_text(size = 9), 
            axis.text = element_text(size = 9), 
            axis.title = element_text(size = 10), 
            plot.margin = unit(c(8.5,8.5,8.5,8.5), "pt")) + 
      scale_y_continuous(labels = scales::percent) 
  }
  
  graph.4 <- graph.4 + ggplot(sub.graph.data.4, aes(fill=Templated, y=count, x=Position)) + 
    ggtitle(paste(group.name, "(unique tag)", sep=" ")) +
    geom_bar(position="stack", stat="identity") +
    xlab("Position") +
    ylab("Tag")+
    theme_bw() + 
    theme(plot.title = element_text(size = 12.5, face = "bold"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 9), 
          axis.text = element_text(size = 9), 
          axis.title = element_text(size = 10), 
          plot.margin = unit(c(8.5,8.5,8.5,8.5), "pt")) 
}
graph.4 <- graph.4 + 
  plot_layout(ncol = 2, guides = "collect")  
ggsave(file=paste(args[[1]], 'graph_4.pdf', sep=''), plot = graph.4, width = 15 * type_ext.scaling, height = 10 * group.scaling, limitsize = FALSE)


# Visualise graph 5
graph.data.5 <- read.csv(paste(args[[2]], 'graph_5_data.csv', sep=''))
names(graph.data.5)[names(graph.data.5) == "nucleotide"] <- "Nucleotide"
names(graph.data.5)[names(graph.data.5) == "group"] <- "Group"
names(graph.data.5)[names(graph.data.5) == "position"] <- "Position"

grouped.graph.data.5 <- graph.data.5 %>% group_by(Group) %>% group_data() 
grouped.graph.data.5 <- as.data.frame(grouped.graph.data.5)
graph.5 <- NULL
group.scaling <- length(unique(graph.data.5[, 'Group'])) / 2 
type_ext.scaling <- length(unique(graph.data.5[, 'Position'])) / 10 

for (row in 1:nrow(grouped.graph.data.5)) {
  group.name = grouped.graph.data.5[row,'Group']
  indexes = unlist(grouped.graph.data.5[row,'.rows'])
  sub.graph.data.5 <- graph.data.5[indexes,]
  
  sub.graph.data.5$Position=factor(sub.graph.data.5$Position, levels = unique(sub.graph.data.5$Position))
  
  if (is.null(graph.5)) {
    graph.5 <- ggplot(sub.graph.data.5, aes(fill=Nucleotide, y=count, x=Position)) + 
      ggtitle(paste(group.name, "(%)", sep=" ")) +
      geom_bar(position="fill", stat="identity")+
      xlab("Position") +
      ylab("%") +
      theme_bw() + 
      theme(
        plot.title = element_text(size = 12.5, face = "bold"), 
        legend.title = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 9), 
        axis.text = element_text(size = 9), 
        axis.title = element_text(size = 10),
        plot.margin = unit(c(8.5,8.5,8.5,8.5), "pt")) + 
      scale_y_continuous(labels = scales::percent) 
  } else {
    graph.5 <- graph.5 + ggplot(sub.graph.data.5, aes(fill=Nucleotide, y=count, x=Position)) + 
      ggtitle(paste(group.name, "(%)", sep=" ")) +
      geom_bar(position="fill", stat="identity") +
      xlab("Position") +
      ylab("%") +
      theme_bw() + 
      theme(
        plot.title = element_text(size = 12.5, face = "bold"), 
        legend.title = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 9), 
        axis.text = element_text(size = 9), 
        axis.title = element_text(size = 10), 
        plot.margin = unit(c(8.5,8.5,8.5,8.5), "pt")) + 
      scale_y_continuous(labels = scales::percent) 
  }
  
  graph.5 <- graph.5 + ggplot(sub.graph.data.5, aes(fill=Nucleotide, y=count, x=Position)) + 
    ggtitle(paste(group.name, "(unique tag)", sep=" ")) +
    geom_bar(position="stack", stat="identity") +
    xlab("Position") +
    ylab("Tag") +
    theme_bw() + 
    theme(plot.title = element_text(size = 12.5, face = "bold"), 
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 9), 
          axis.text = element_text(size = 9), 
          axis.title = element_text(size = 10), 
          plot.margin = unit(c(8.5,8.5,8.5,8.5), "pt")) 
}
graph.5 <- graph.5 + 
  plot_layout(ncol = 2, guides = "collect") 
ggsave(file=paste(args[[1]], 'graph_5.pdf', sep=''), plot = graph.5, width = 10 * type_ext.scaling, height = 10 * group.scaling, limitsize = FALSE)

# Visualise graph 6
graph.data.6 <- read.csv(paste(args[[2]], 'graph_6_data.csv', sep=''))

names(graph.data.6)[names(graph.data.6) == "templated"] <- "Templated"
names(graph.data.6)[names(graph.data.6) == "group"] <- "Group"
names(graph.data.6)[names(graph.data.6) == "position"] <- "Position"
names(graph.data.6)[names(graph.data.6) == "type"] <- "Type"

grouped.graph.data.6 <- graph.data.6 %>% group_by(Group,Type) %>% group_data() 
grouped.graph.data.6 <- as.data.frame(grouped.graph.data.6)
graph.6 <- NULL

group.scaling <- length(unique(graph.data.6[, 'Group'])) / 2 
type.scaling <- length(unique(graph.data.6[, 'Position'])) / 10 

for (row in 1:nrow(grouped.graph.data.6)) {
  group.name = grouped.graph.data.6[row,'Group']
  type.name = grouped.graph.data.6[row,'Type']
  
  indexes = unlist(grouped.graph.data.6[row,'.rows'])
  sub.graph.data.6 <- graph.data.6[indexes,]
  
  sub.graph.data.6 <- sub.graph.data.6 %>%
    mutate(int.position = as.integer(Position)) 
  
  max_position <- sub.graph.data.6 %>%
    filter(!is.na(int.position)) %>%
    group_by(int.position) %>%
    summarize(total_count = sum(count)) %>%
    filter(total_count > 0) %>%
    summarize(max_position = max(int.position)) %>%
    pull(max_position)

  sub.graph.data.6 <- sub.graph.data.6 %>% 
    filter(is.na(int.position) | int.position <= max_position)

  sub.graph.data.6$Position=factor(sub.graph.data.6$Position, levels = unique(sub.graph.data.6$Position))
  
  if (is.null(graph.6)) {
    graph.6 <- ggplot(sub.graph.data.6, aes(fill=Templated, y=count, x=Position)) + 
      ggtitle(paste(group.name, type.name, "(unique tag)", sep=" ")) +
      geom_bar(position="stack", stat="identity") +
      xlab("Position") +
      ylab("Tag")+
      theme_bw() + 
      theme(plot.title = element_text(size = 12.5, face = "bold"), 
            legend.title = element_text(size = 10, face = "bold"), 
            legend.text = element_text(size = 9), 
            axis.text = element_text(size = 9), 
            axis.title = element_text(size = 10), 
            plot.margin = unit(c(8.5,8.5,8.5,8.5), "pt")) 
  } else {
    graph.6 <- graph.6 + ggplot(sub.graph.data.6, aes(fill=Templated, y=count, x=Position)) + 
      ggtitle(paste(group.name, type.name, "(unique tag)", sep=" ")) +
      geom_bar(position="stack", stat="identity") +
      xlab("Position") +
      ylab("Tag")+
      theme_bw() + 
      theme(plot.title = element_text(size = 12.5, face = "bold"), 
            legend.title = element_text(size = 10, face = "bold"), 
            legend.text = element_text(size = 9), 
            axis.text = element_text(size = 9), 
            axis.title = element_text(size = 10), 
            plot.margin = unit(c(8.5,8.5,8.5,8.5), "pt")) 
  }
}  
graph.6 <- graph.6 + 
  plot_layout(ncol = 2, guides = "collect") 
ggsave(file=paste(args[[1]], 'graph_6.pdf', sep=''), plot = graph.6, width = 10 * type.scaling, height = 10 * group.scaling, limitsize = FALSE)



















