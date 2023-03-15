library(ggplot2)
library(gtable)
library(grid)


hypo<-data.frame(Genotype=c("N","O","N","O"), Site=c("O","O","N","N"),
                 Response=c(7,8.3,8,4))

hypo$Genotype<-factor(hypo$Genotype, levels=c("O","N"))
hypo$Site<-factor(hypo$Site, levels=c("O","N"))
ggplot(hypo, aes(x=Site, y=Response, fill= Genotype, color=Genotype))+
    geom_bar(position=position_dodge(.6), stat="identity",width=0.5)+
    ylab('Reproductive success')+
    scale_x_discrete(labels=c("Clean water", "Poor water quality"))+
    #scale_y_continuous(breaks='none')+
    scale_fill_manual (values=c("#4E74FF","#F4CC5C"))+
    scale_color_manual(values=c("#4E74FF","#F4CC5C"))+
    theme_bw()+xlab('')+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          legend.title = element_blank(),axis.text.x = element_text(size=11, color=1),
          legend.text = element_text(size=12, color=1),
          axis.ticks=element_blank(), axis.text.y = element_blank(),
          panel.border = element_blank(),
          axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last")),
          axis.line.x = element_line())
ggsave("Output/EnvxGeno_hypo_plot.pdf", width = 4, height = 3)
          