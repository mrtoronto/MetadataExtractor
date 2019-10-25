# Distribution of ages from GEO sample metadata
# ... note this figure highlights errors in age extraction, exclusions are made for now
#
# Author: Seth Rhoades

library(jsonlite, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(ggthemes, quietly = TRUE)
library(readr, quietly = TRUE)

filePath = file.path('.')
refDirectory = c('refFiles')
metaFile = c('GEO_MusmusculusMetadata.json')

metaData = fromJSON(txt = readLines(file.path(filePath, refDirectory, metaFile)))
ages = c()
for(sample in metaData){
    if(sample$Age != c('n/a')){
        if(sample$Age < 125){
            ages = c(ages, sample$Age)
        } else {print(paste('Suspected age error', sample$ID, ':', sample$Age))}     
    }
}
ages = as.data.frame(ages)
colnames(ages) = c('Age')

ageHist = ggplot(ages, aes(x = Age)) + 
        geom_histogram(colour = '#303030', binwidth = 4) + 
        theme_light() +
        ylab('Count') + 
        xlab('Age (weeks)')+
        ggtitle('Mus musculus age estimates: NIH gene expression samples') + 
        theme(legend.position = 'none', 
            axis.text = element_text(size = 10),
            title = element_text(size = 11, face = "bold"), 
            axis.title = element_text(size = 11, face = "bold"),
            panel.grid = element_line(colour = 'grey90'),
            plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill = '#FFFDE9'),
            panel.grid.major.x = element_line(colour = 'grey92', size = 0.4),
            panel.grid.major.y = element_line(colour = 'grey92', size = 0.4),
            panel.grid.minor.x = element_line(colour = 'grey92', size = 0.4),
            panel.grid.minor.y = element_line(colour = 'grey92', size = 0.4))

ggsave('Figures/MusMusculusAgeHistogram.png', width = 6, height = 5, dpi = 600)
