library(ggplot2)
library(stringr)

banff_score_table <- read.table('~/Desktop/histomx/test_files//banff_score_table.txt', sep='\t', header=TRUE, check.names=FALSE)
rownames(banff_score_table) <- banff_score_table$score

##-----------------
## stacked bar plot

df.banff <- data.table::melt(banff_score_table)
df.banff$variable <- factor(df.banff$variable, levels=c("3", "2", "1", "0"))
#df.banff$score <- factor(df.banff$score, levels=c("g", "ptc", "i", "t"))
df.banff$score <- factor(df.banff$score, levels=c("v", "t", "i", "ptc", "g"))

#my_cols=c("#EFF3FF", "#BDD7E7", "#6BAED6", "darkblue")
my_cols=c("steelblue4", "steelblue", "lightblue", "lightcyan")
#my_cols=c("steelblue4", "steelblue", "#BDD7E7", "#EFF3FF")
	
ggplot(data=df.banff) + ylab("probability (%)") + xlab("Banff score") +
	geom_bar(aes(x=score, y=value, fill=variable), color="black", stat="identity", alpha=1, width=0.7) + 
	scale_fill_manual("score", values=my_cols) +
	coord_flip() +
	theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
	      panel.background=element_blank(), axis.line=element_line(colour="white"), 
	      axis.text=element_text(size=20, color="navy"), 
	      axis.title.x=element_text(size=18, color="darkslateblue"),
	      axis.title.y=element_blank(),
	      legend.title=element_text(size=14, color="navy"),
	      legend.text=element_text(size=12, color="navy"), legend.key.size = unit(1, 'cm'),
	      panel.border=element_rect(colour="gray", fill=NA, size=1))

##-----------------
## radar chart
radar_table <- banff_score_table
radar_table$score <- NULL
tab <- data.frame(max=apply(radar_table, 1, max),
		  max_score=colnames(radar_table)[max.col(radar_table)])
tab$max_score <- paste0(rownames(tab), tab$max_score)

ggplot(tab) +
	geom_hline(aes(yintercept = y), data.frame(y = c(0:100) ), color = "white") + 
	geom_col(aes(x = reorder(str_wrap(max_score, 7), max), y = max, fill = max),
		 position = "dodge2", show.legend = TRUE, alpha = .9) +
	geom_segment(aes(x = reorder(str_wrap(max_score, 7), max), y = 0,
			 xend = reorder(str_wrap(max_score, 7), max),
			 yend = max(max)), linetype = "dashed", color = "navy") + 
	coord_polar(clip="off") +
	scale_fill_gradientn(colours=rev(my_cols), limits=c(0,100)) +
	theme(axis.title = element_blank(),
	      axis.ticks = element_blank(),
	      axis.text.y = element_blank(),
	      axis.text.x = element_text(color="navy", size=24, face="bold"),
	      panel.background = element_rect(fill = 'white', colour = 'white'),
	      #plot.background=element_rect(fill="white"),
	      #panel.grid.minor=element_line(colour="gray"),
	      panel.grid.major=element_line(colour="gray"),
	      panel.border=element_rect(colour=NA, fill=NA, size=5),
	      legend.title = element_blank(),
	      legend.text = element_text(size=14, color="navy"), legend.key.size=unit(1.5, 'cm'),
	      legend.position="right")

##-----------------
## bubble chart
df.banff$variable <- factor(df.banff$variable, levels=c("0", "1", "2", "3"))

ggplot(df.banff, aes(x=variable, y=score, size=value)) + xlab("") + ylab("") +
	geom_point(aes(fill=value), pch=21, colour="azure4", alpha=1) +
	scale_size(range = c(0, 30), name="%", guide="none") + #increase relative bubble size
	scale_fill_gradientn(colours=rev(my_cols), limits=c(0,100)) +
	theme(panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor=element_blank(),
	      panel.background=element_blank(), axis.line=element_line(colour="white"),
	      axis.text=element_text(size=20, color="navy"),
	      axis.title.x=element_text(size=20, color="navy"), 
	      axis.title.y=element_text(size=20, color="navy"),
	      legend.title=element_blank(),
	      legend.text=element_text(size=16, color="darkslateblue"), legend.key.size = unit(1.5, 'cm'),
	      panel.border=element_rect(colour="gray", fill=NA, size=1))

##-----------------
## heatmap
df.banff$value <- round(df.banff$value, 1)

ggplot(df.banff, aes(x=variable, y=score, fill=value)) + xlab("") + ylab("") +
	geom_tile() + scale_fill_gradientn(colours=rev(my_cols), limits=c(0,100)) +
	geom_text(aes(label=value), size=5, color="navy") +
	theme(panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor=element_blank(),
	      panel.background=element_blank(), axis.line=element_line(colour="white"),
	      axis.text=element_text(size=20, color="navy"),
	      axis.title.x=element_text(size=20, color="navy"), 
	      axis.title.y=element_text(size=20, color="navy"),
	      legend.title=element_blank(),
	      legend.text=element_text(size=16, color="darkslateblue"), legend.key.size = unit(1.5, 'cm'),
	      panel.border=element_rect(colour="gray", fill=NA, size=1))


