library(ggplot2)
library(scales)
library(dplyr)
library(gridExtra)

# create graph folder 
dir.create(file.path('../data/graphs/'))

# Visualise graph 1
sub.graphs.1 <- list()
graph.data.1 <- read.csv('../data/8_graph_processed_data/graph_1_data.csv')
graph.data.1$group<-as.factor(graph.data.1$group)
sub.graphs.1[[1]] <- ggplot(data = graph.data.1, aes(x = group, y = rpm, group = type, color = type)) +
  geom_line() + 
  geom_point() + 
  theme_bw()+
  theme(legend.text = element_text(size=14),
        legend.title = element_blank(),
        axis.title=element_text(size=14), 
        axis.text = element_text(size=14))
sub.graphs.1[[2]] <- ggplot(data = graph.data.1, aes(x = group, y = relative_abundance, group = type, color = type)) +
  geom_line() + 
  geom_point() + 
  theme_bw()+
  theme(legend.text = element_text(size=14),
        legend.title = element_blank(),
        axis.title=element_text(size=14), 
        axis.text = element_text(size=14))
#do.call('grid.arrange', c(sub.graphs.1, ncol=2))
graph.1 <- grid.arrange(grobs=sub.graphs.1, ncol = 2)
ggsave(file='../data/graphs/graph_1.pdf', graph.1)

# Visualise graph 2
sub.graphs.2 <- list()
graph.data.2 <- read.csv('../data/8_graph_processed_data/graph_2_data.csv')
grouped.graph.data.2 <- graph.data.2 %>% group_by(group) %>% group_data() 
grouped.graph.data.2 <- as.data.frame(grouped.graph.data.2)
graph.count.2 = 1
for (row in 1:nrow(grouped.graph.data.2)) {
  group.name = grouped.graph.data.2[row,'group']
  indexes = unlist(grouped.graph.data.2[row,'.rows'])
  sub.graph.data.2 <- graph.data.2[indexes,]
  
  sub.graphs.2[[graph.count.2]] <- ggplot(sub.graph.data.2, aes(x="", y=rpm, fill=grouped_type)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() +
    geom_text(aes(label = paste(round(rpm/sum(rpm)*100,1), "%")),
              position = position_stack(vjust = 0.5), color="white")+
    ggtitle(paste(group.name, "isomiRNAs (rpm)", sep=" "))+
    theme(plot.title = element_text(hjust = 0.5))
  graph.count.2 = graph.count.2 + 1 
  
  sub.graphs.2[[graph.count.2]] <- ggplot(sub.graph.data.2, aes(x="", y=unique_tag, fill=grouped_type)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() +
    geom_text(aes(label = paste(round(unique_tag/sum(unique_tag)*100,1), "%")),
              position = position_stack(vjust = 0.5), color="white")+
    ggtitle(paste(group.name, "isomiRNAs (unique tag)", sep=" "))+
    theme(plot.title = element_text(hjust = 0.5))
  graph.count.2 = graph.count.2 + 1 
}
graph.2 <- grid.arrange(grobs=sub.graphs.2, ncol = 2)
ggsave(file='../data/graphs/graph_2.pdf', graph.2)

# Visualise graph 3
sub.graphs.3 <- list()
graph.data.3 <- read.csv('../data/8_graph_processed_data/graph_3_data.csv')
grouped.graph.data.3 <- graph.data.3 %>% group_by(group) %>% group_data() 
grouped.graph.data.3 <- as.data.frame(grouped.graph.data.3)
graph.count.3 = 1
for (row in 1:nrow(grouped.graph.data.3)) {
  group.name = grouped.graph.data.3[row,'group']
  indexes = unlist(grouped.graph.data.3[row,'.rows'])
  sub.graph.data.3 <- graph.data.3[indexes,]
  
  sub.graphs.3[[graph.count.3]] <- ggplot(data=sub.graph.data.3, aes(x=group, y=rpm, fill=type_nt)) + 
    geom_bar(position = "fill", stat="identity") + 
    theme_bw()+
    theme(axis.ticks.x=element_blank()) + 
    theme(panel.grid.major.x = element_blank())+ 
    theme(panel.grid.major.y = element_blank()) + 
    ylab("") +  
    xlab(" ") +
    theme(axis.title.y = element_text(size = 11, colour = "black"))+ 
    theme(axis.text.y = element_text(size = 11, colour = "black")) +
    theme(plot.title = element_text(hjust=0.5, face = "bold", size=10))+ 
    scale_y_continuous(labels = scales::percent) + 
    theme(axis.text.y = element_text(size = 11, colour = "black"))
  graph.count.3 = graph.count.3 + 1 
  
  sub.graphs.3[[graph.count.3]] <- ggplot(data=sub.graph.data.3, aes(x=group, y=unique_tag, fill=type_nt)) + 
    geom_bar(position = "fill", stat="identity") + 
    theme_bw()+
    theme(axis.ticks.x=element_blank()) + 
    theme(panel.grid.major.x = element_blank())+ 
    theme(panel.grid.major.y = element_blank()) + 
    ylab("") +  
    xlab(" ") +
    theme(axis.title.y = element_text(size = 11, colour = "black"))+ 
    theme(axis.text.y = element_text(size = 11, colour = "black")) +
    theme(plot.title = element_text(hjust=0.5, face = "bold", size=10))+ 
    scale_y_continuous(labels = scales::percent) + 
    theme(axis.text.y = element_text(size = 11, colour = "black"))
  graph.count.3 = graph.count.3 + 1 
}
graph.3 <- grid.arrange(grobs=sub.graphs.3, ncol = 2)
ggsave(file='../data/graphs/graph_3.pdf', graph.3)

# Visualise graph 4
sub.graphs.4 <- list()
graph.data.4 <- read.csv('../data/8_graph_processed_data/graph_4_data.csv')
grouped.graph.data.4 <- graph.data.4 %>% group_by(group) %>% group_data() 
grouped.graph.data.4 <- as.data.frame(grouped.graph.data.4)
graph.count.4 = 1
for (row in 1:nrow(grouped.graph.data.4)) {
  group.name = grouped.graph.data.4[row,'group']
  indexes = unlist(grouped.graph.data.4[row,'.rows'])
  sub.graph.data.4 <- graph.data.4[indexes,]
  
  sub.graph.data.4$Position=factor(sub.graph.data.4$Position, levels = unique(sub.graph.data.4$Position))
  
  sub.graphs.4[[graph.count.4]] <- ggplot(sub.graph.data.4, aes(fill=Templated, y=count, x=Position)) + 
    geom_bar(position="fill", stat="identity")+
    xlab("Position") +
    ylab("%")+
    theme_bw()+
    theme(legend.text = element_text(size=14),
          legend.title = element_text(size=14),
          axis.title=element_text(size=14), 
          axis.text = element_text(size=14))
  graph.count.4 = graph.count.4 + 1 
  
  sub.graphs.4[[graph.count.4]] <- ggplot(sub.graph.data.4, aes(fill=Templated, y=count, x=Position)) + 
    geom_bar(position="stack", stat="identity")+
    xlab("Position") +
    ylab("Tag")+
    theme_bw()+
    theme(legend.text = element_text(size=14),
          legend.title = element_text(size=14),
          axis.title=element_text(size=14), 
          axis.text = element_text(size=14))
  graph.count.4 = graph.count.4 + 1 
}
graph.4 <- grid.arrange(grobs=sub.graphs.4, ncol = 2)
ggsave(file='../data/graphs/graph_4.pdf', graph.4)

# Visualise graph 4
sub.graphs.5 <- list()
graph.data.5 <- read.csv('../data/8_graph_processed_data/graph_5_data.csv')
grouped.graph.data.5 <- graph.data.5 %>% group_by(group) %>% group_data() 
grouped.graph.data.5 <- as.data.frame(grouped.graph.data.5)
graph.count.5 = 1
for (row in 1:nrow(grouped.graph.data.5)) {
  group.name = grouped.graph.data.5[row,'group']
  indexes = unlist(grouped.graph.data.5[row,'.rows'])
  sub.graph.data.5 <- graph.data.5[indexes,]
  
  sub.graph.data.5$Position=factor(sub.graph.data.5$Position, levels = unique(sub.graph.data.5$Position))
  
  sub.graphs.5[[graph.count.5]] <- ggplot(sub.graph.data.5, aes(fill=Nucleotide, y=count, x=Position)) + 
    geom_bar(position="fill", stat="identity")+
    xlab("Position") +
    ylab("%")+
    theme_bw()+
    theme(legend.text = element_text(size=14),
          legend.title = element_text(size=14),
          axis.title=element_text(size=14), 
          axis.text = element_text(size=14))
  graph.count.5 = graph.count.5 + 1 
  
  sub.graphs.5[[graph.count.5]] <- ggplot(sub.graph.data.5, aes(fill=Nucleotide, y=count, x=Position)) + 
    geom_bar(position="stack", stat="identity")+
    xlab("Position") +
    ylab("Tag")+
    theme_bw()+
    theme(legend.text = element_text(size=14),
          legend.title = element_text(size=14),
          axis.title=element_text(size=14), 
          axis.text = element_text(size=14))
  graph.count.5 = graph.count.5 + 1 
}
graph.5 <- grid.arrange(grobs=sub.graphs.5, ncol = 2)
ggsave(file='../data/graphs/graph_5.pdf', graph.5)

