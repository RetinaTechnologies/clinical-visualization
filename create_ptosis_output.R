# Install packages if not exist
package_list=c('data.table', 'gridExtra', 'ggplot2', 'cowplot', 'plyr', 
				 'dplyr', 'grid', 'pracma')
new_packages=package_list[!(package_list %in% installed.packages()[,'Package'])]
if(length(new_packages)) install.packages(new_packages)

# Load packages
library(data.table)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(plyr)
library(dplyr)
library(grid)
library(pracma)

# TODO: Set path to code
setwd('')

# Define figure theme
theme_basic <-  function(axis_size=0.5, title_size=8, subtitle_size=6,
                  col_gen='grey50', legend_title_size=0.5, 
                  legend_text_size=0.4, legend_tick_size=0.08,
                  legend_width=0.5, legend_height=0.2, 
                  legend_hjust_title=0.5)  {

  theme_bw() +
  theme(
    plot.title=element_text(size=title_size, colour=col_gen, face='bold'),
    plot.subtitle=element_text(size=subtitle_size, colour=col_gen,
      face='plain'),
    plot.caption=element_text(size=(subtitle_size-1), colour=col_gen,
      face='plain'),
    legend.position='bottom',
    legend.key.height=unit(legend_height,'cm'),
    legend.key.width=unit(legend_width,'cm'),
    axis.ticks.length=unit(legend_tick_size,'cm'),
    legend.title=element_text(size=rel(legend_title_size), colour=col_gen,
      hjust=legend_hjust_title, face='plain'),
    legend.text=element_text(size=rel(legend_text_size), colour=col_gen), 
    panel.grid.major=element_blank(), 
    panel.grid.minor=element_blank(), 
    axis.title.x=element_blank(), 
    axis.text.x=element_blank(), 
    axis.ticks.x=element_blank(), 
    axis.title.y=element_blank(), 
    axis.text.y=element_blank(), 
    axis.ticks.y=element_blank()    
   )
}

# Calculate area under the curve
auc <- function(y){
  n <- length(y)
  0.5*(y[1]+y[n]+2*sum(y[-c(1,n)]))
}

create_ptosis_output <- function(patient, eye, untaped_only=FALSE) {
	# Set path for files containing taped and untaped exams
	if(!untaped_only) {
		file_taped <- paste0('../data/modified/', patient, 
							 '_Superior_', eye, '_taped_eye_mod.csv')
	}

	file_untaped <- paste0('../data/modified/', patient, 
						   '_Superior_', eye, '_eye_mod.csv')
	if(untaped_only) {
		lst <- list(file_untaped)
	} else {
		lst <- list(file_untaped, file_taped)
	}

	# Generate random patient identifier
	patient_id <- paste(sample(0:9, 7, replace=TRUE), collapse="")
	pdf(paste0('../output/', patient_id, '_Superior_', eye, '_eye.pdf'))
	for (file in lst) {
		print(file)
		dt <- fread(file)
		
		# Preprocess data table
		names(dt) <- make.names(names(dt), unique=TRUE)

		# Filter out untested points
		dt_tested <- dt[Was.point.tested == TRUE]

		# Determine margin
		dt_hits <- dt_tested[Response.time.s. > 0]
		dt_hits_max <- data.table(dt_hits %>% group_by(x) %>% top_n(n=1, wt=y))
		dt_misses <- setdiff(dt_tested, dt_hits)

		# Plot scatter without best fit
		sctr_no_line <- ggplot() + 
			    		geom_point(dt_hits, mapping=aes(x, y), shape=8, size=2) +
			    		geom_point(dt_misses, mapping=aes(x, y), shape=15, size=2) +
			    		geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
			    		geom_hline(yintercept=3, linetype='dashed') +
			    		xlim(-3, 3) + ylim(-4, 4) + 
			    		annotate('text', x=-3, y=4, label=paste0('Patient: ', patient_id), hjust=0) +
			    		annotate('text', x=-3, y=3.8, label=paste0('Eye: ', eye), hjust=0) + 
			    		theme_basic()		
		# Plot scatter with best fit
		sctr <- ggplot() + 
			    geom_point(dt_hits, mapping=aes(x, y), shape=8, size=2) +
			    geom_point(dt_misses, mapping=aes(x, y), shape=15, size=2) +
			    stat_smooth(dt_hits_max, mapping=aes(x, y), method='auto', se=FALSE, 
			    			color='black', size=0.8, linetype=1) +
			    geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
			    geom_hline(yintercept=3, linetype='dashed') +
			    xlim(-3, 3) + ylim(-4, 4) + 
			    annotate('text', x=-3, y=4, label=paste0('Patient: ', patient_id), hjust=0) +
			    annotate('text', x=-3, y=3.8, label=paste0('Eye: ', eye), hjust=0) + 
			    theme_basic()
		auc_top=round(trapz(rep(3, 11)), 2)
		
		# Include annotations and calculate area under the curve
		if(grepl('taped', file)) {
			auc_taped <- round(trapz(dt_hits_max$y), 2)
			percent_change <- round((auc_taped - auc_untaped)/auc_taped * 100, 2)
			sctr_no_line <- sctr_no_line + annotate('text', x=-3, y=3.6, label='Taped: Yes', hjust=0) + 
			annotate('text', x=-3, y=3.4, label=paste0('AUC: ', auc_taped), hjust=0) + 
			annotate('text', x=-3, y=3.2, label=paste0('Percent change: ', percent_change, '%'), hjust=0)			
			sctr <- sctr + annotate('text', x=-3, y=3.6, label='Taped: Yes', hjust=0) + 
			annotate('text', x=-3, y=3.4, label=paste0('AUC: ', auc_taped), hjust=0) + 
			annotate('text', x=-3, y=3.2, label=paste0('Percent change: ', percent_change, '%'), hjust=0)
		} else {
			auc_untaped=round(trapz(dt_hits_max$y), 2)
			sctr_no_line <- sctr_no_line + annotate('text', x=-3, y=3.6, label='Taped: No', hjust=0) + 
			annotate('text', x=-3, y=3.4, label=paste0('AUC: ', auc_untaped), hjust=0) +
			annotate('text', x=-3, y=3.2, label=paste0('AUC max: ', auc_top), hjust=0)			
			sctr <- sctr + annotate('text', x=-3, y=3.6, label='Taped: No', hjust=0) + 
			annotate('text', x=-3, y=3.4, label=paste0('AUC: ', auc_untaped), hjust=0) +
			annotate('text', x=-3, y=3.2, label=paste0('AUC max: ', auc_top), hjust=0)
		}
		print(sctr_no_line)
		print(sctr)
	}
	dev.off()
}

# Command line args specifying patient and eye for which to create output
args=commandArgs(trailingOnly=TRUE)
if(length(args) == 0) {
	stop('Must specify [patient] and [eye]')
} else {
	patient=args[1]
	eye=args[2]
}
create_ptosis_output(patient, eye)