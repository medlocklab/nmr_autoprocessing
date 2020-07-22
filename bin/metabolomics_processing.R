library(plyr)
library(ggplot2)
library(pheatmap)
library(multcomp)
#

# Load integrals and sample data and remove metabolites we don't need ranges for
path <- "~/Documents/projects/asf_modeling/data/metabolomics/coculture_metabolomics/"
integral_file <- 'integrals_20170323.txt'
integral_range_file <- 'coculture_peak_ranges.txt'
integrals = read.table(paste0(path,integral_file),header=FALSE,
                         sep = ",")
integral_ranges = read.table(paste0(path,integral_range_file),header=TRUE,
                             sep = "\t")

master <- read.table(file=paste0(path,"merged_metadata.txt"),header=TRUE,sep='\t')
#


# set column names of integrals using integral_ranges met names
colnames(integrals) <- integral_ranges$met_name
mets = colnames(integrals)
# Set rownames using experiment_id. This should be ordered correctly.
integrals$experiment_id <- master$experiment_id
# grab metname list for calculating z-scores


# merge so that classes and integrals are in one dataframe
all_data <- join(master,integrals,by="experiment_id")

# get mean and std dev for each metabolite from blanks
blanks = subset(all_data,species == '0')
blank_means = colMeans(blanks[,names(blanks) %in% mets])
blank_stdev = apply(blanks[,names(blanks) %in% mets],2,sd)

subtracted_integrals <- sweep(integrals[,!names(integrals)=='experiment_id'], 2, blank_means, `-`)
z_scores <- sweep(subtracted_integrals, 2, blank_stdev, `/`)

# replace integral values with z-scores
all_data[,mets] <- z_scores

# Find negative values and replace them with the lowest non-negative value for that metabolite
log_integrals = integrals
log_integrals[,colnames(log_integrals) %in% mets] = 
        (apply(log_integrals[,colnames(log_integrals) %in% mets],2,function(x){x[x<0] = min(x[x>0]);x}))
# generate z-scores using log-transformed data
log_all_data <- join(master,log_integrals,by="experiment_id")
log_blanks = subset(log_all_data,species == '0')
log_integrals[,mets] = log(log_integrals[,mets])
log_blank_means = colMeans(log_blanks[,names(log_blanks) %in% mets])
log_blank_stdev = apply(log_blanks[,names(log_blanks) %in% mets],2,sd)
log_subtracted_integrals <- sweep(log_integrals[,!names(log_integrals)=='experiment_id'], 2, log_blank_means, `-`)
log_z_scores <- sweep(log_subtracted_integrals, 2, log_blank_stdev, `/`)

# replace integral values with z-scores
log_all_data[,mets] <- log_z_scores

# remove experiment 4
log_all_data = log_all_data[log_all_data$Run != "4",]

# Test for differential abundance
species_list = c('356','360','361','492','500','502','519')
for (species_tested in species_list) {
  # subset the master dataframe by cultures that included species
  contains_species = log_all_data[grepl(species_tested,log_all_data$species),]
  contains_species = droplevels(contains_species)
  # get the correct qPCR probe for the species being tested
  qPCR_column = colnames(contains_species)[grepl(species_tested,colnames(contains_species))]
  unique_combos = as.character(unique(contains_species$species)[!(unique(contains_species$species) %in% c(species_tested))])
  # for each unique combo, perform the t test
  x = contains_species[contains_species$species == species_tested,][,qPCR_column]
  results = lapply(seq_along(unique_combos), function (n) {
    y = contains_species[contains_species$species == unique_combos[n],][,qPCR_column]
    result = t.test(x,y)
    return(result)})
  names(results) <- paste(matrix(species_tested, ncol = 2, byrow = TRUE)[,1], matrix(unlist(unique_combos), ncol = 1, byrow = TRUE), sep = " vs. ")
  # generate vector of p values for multiple testing correction
  p_init = numeric()
  p_vector = c(p_init)
  t_stat_init = numeric()
  t_stat_vector = c(t_stat_init)
  for (i in 1:length(names(results))) {
    p_vector = c(p_vector,results[[names(results)[i]]]$p.value)
    t_stat_vector = c(t_stat_vector,results[[names(results)[i]]]$statistic)
  }
  names(p_vector) = names(results)
  names(t_stat_vector) = names(results)
  # perform the correction
  final_p = p.adjust(p_vector,method="BH")
  print(t_stat_vector)
  print(final_p)
  resultframe = data.frame(final_p,t_stat_vector)
  # save the results
  write.table(resultframe,file=paste0(path,paste0(species_tested,"_diff_abundance.txt")),row.names=TRUE,sep='\t')
}

# calculate metabolite yields for each mono-culture sample by dividing the z-score by species abundance.
species_list = c('356','360','361','492','500','502','519')

#first remove experiment 4
all_data = all_data[all_data$Run != "4",]
single_species_all = all_data[all_data$species %in% species_list,]
qpcr_cols = colnames(single_species_all)[grep('qpcr',colnames(single_species_all))]
# match the species names to the qpcr column
#single_species_all$qpcr_name = NA
for (species_name in species_list) {
  #single_species_all[single_species_all$species == species_name,]$qpcr_name = qpcr_cols[grep(species_name,qpcr_cols)]
  qpcr_name = qpcr_cols[grep(species_name,qpcr_cols)]
  single_species_all[single_species_all$species == species_name,][,mets] = single_species_all[single_species_all$species == species_name,][,mets]/single_species_all[single_species_all$species == species_name,][,qpcr_name]
}

avg_single_yield = aggregate(single_species_all[,mets], list(single_species_all$species), mean)
rownames(avg_single_yield) = avg_single_yield$Group.1
avg_single_yield$Group.1 = NULL

copy_all_data = all_data
# replace NA with 0
copy_all_data[is.na(copy_all_data)] = 0
copy_all_data[,mets] = 0
for (species_name in species_list) {
  qpcr_name = qpcr_cols[grep(species_name,qpcr_cols)]
  copy_all_data[,mets] = copy_all_data[,mets] + as.matrix(copy_all_data[,qpcr_name])%*%as.matrix(avg_single_yield[species_name,])
}

# for each co-culture, generate a predicted metabolite abundance using the following, where
# all yields are the average across all monoculture samples:
# predicted = yield[species1]*abundance[species1] + yield[species2]*abundance[species2]

# now compare each co-culture group in copy_all_data (predictions) to all_data (actual values)
# for each metabolite, and report the p adj and t stat.
all_combos = unique(copy_all_data$species)
all_combos = all_combos[!(all_combos %in% species_list)]
all_combos = all_combos[!(all_combos %in% c('0'))]
all_combos = droplevels(all_combos)

#only test known metabolites
#known_mets = mets[-grep("unknown", mets)]


all_met_tests = NULL
for (combo in all_combos) {
  subset_predicted = copy_all_data[copy_all_data$species == combo,mets]
  subset_observed = all_data[all_data$species == combo,mets]
  combo_p_init = numeric()
  combo_ps = (combo_p_init)
  combo_t_init = numeric()
  combo_ts = c(combo_t_init)
  met_order_init = character()
  met_order = c(met_order_init)
  for (met in mets) {
    result = t.test(subset_observed[,met],subset_predicted[,met])
    #f (result$p.value < 0.01) {
    #  print(combo)
    #  print(met)
    #  print(result$p.value)
    #}
    combo_ps = c(combo_ps,result$p.value)
    combo_ts = c(combo_ts,result$statistic)
    met_order = c(met_order,met)
  }
  
  # make new df from p,t,met, and combo
  collapse = data.frame(combo_ps,combo_ts,met_order,combo)
  all_met_tests = rbind(collapse,all_met_tests)
  
}

# FDR correct the p values
all_met_tests$combo_ps = p.adjust(all_met_tests$combo_ps,method="BH")
# save the result
write.table(all_met_tests,file=paste0(path,"all_diff_metabolites.txt"),row.names=FALSE,sep='\t')

# make heatmap for 356 519 comparison


interaction_df = colMeans(all_data[all_data$species == '356',known_mets])
single_519_profile = colMeans(all_data[all_data$species == '519',known_mets])
interaction_df = rbind(interaction_df,single_519_profile)
predicted = colMeans(copy_all_data[copy_all_data$species == '356,519',known_mets])
#interaction_df = rbind(interaction_df,predicted)
observed = colMeans(all_data[all_data$species == '356,519',known_mets])
#interaction_df = rbind(interaction_df,observed)
diff = observed - predicted
interaction_df = rbind(interaction_df, diff)

rownames(interaction_df)[1] = "356 monoculture"
rownames(interaction_df)[2] = "519 monoculture"
#rownames(interaction_df)[1] = "Predicted co-culture"
#rownames(interaction_df)[1] = "Observed co-culture"
rownames(interaction_df)[3] = "Observed - predicted"
  
interaction_df[interaction_df > 10] = 10
interaction_df[interaction_df < -10] = -10
interaction_df[interaction_df < 2 & interaction_df > -2] = 0
pheatmap(interaction_df,filename='~/interaction_356_519.jpg',width=9,height=3.2,cluster_rows=FALSE)

# get sig p vals
subset_356_519 = all_met_tests[all_met_tests$combo == '356,519',]
sig_only = subset_356_519[subset_356_519$combo_ps < 0.05,]$met_order

# some calculations for visualizing the mean vs. std dev distribution.
# compare untransformed to log-transformed
a = log(integrals[,names(integrals) %in% mets])
b = integrals[,names(integrals) %in% mets]
a[is.na(a)] = 0
log_st_dev = apply(a,2,sd)
log_avg = colMeans(a)
normal_st_dev =apply(b,2,sd)
normal_avg = colMeans(b)
log_df = data.frame(log_avg,log_st_dev)
normal_df = data.frame(normal_avg,normal_st_dev)
log_ordered = log_df[order(log_avg),]
log_ordered$rank = seq(nrow(log_ordered))
normal_ordered = normal_df[order(normal_avg),]
normal_ordered$rank = seq(nrow(normal_ordered))

plot(log_ordered$rank,log_ordered$log_st_dev)
plot(normal_ordered$rank,normal_ordered$normal_st_dev)

# generate plots using log-transformed z-scores






#targets = c('500','519','500,519','0')
targets = c('361','519','361,519','0')
targets = c('356','361','519','356,361','356,519','361,519','356,361,519','0')

p <- ggplot(subset(all_data,species %in% targets)) +
        geom_point(aes(species,Lactate),group='species')
p

targets = c('356','360','361','492','500','502','519','0')
single_species = c('519','361','361,519')
single_species = c('356','360','356,360')
single_species = c('356','492','356,492')
single_species = c('356','361','519','356,361','356,519','356,361,519')
single_species = c('360','519','360,519')
single_species = c('356','502','356,502')
single_species = c('356','519','356,519')
single_species = c('500','502','500,502')
single_species = c('356','360','361','492','500','502','519')
#single_species = c('519','0')


# extract only the samples with growth
only_grows = log_all_data[log_all_data$normalized_OD > 0.01,]
# remove an outlier (low sample volume?)
only_grows = only_grows[!only_grows$experiment_id == '122B',]
# collapse as mean
only_grows_means = aggregate(.~species, data=only_grows, mean)

# select species of interest
single_integrals = only_grows_means[only_grows_means$species %in% single_species,]
#non-transformed z-scores
#single_integrals = all_data[all_data$species %in% single_species,]

# remove the unknown metabolites
single_integrals = single_integrals[, -grep("unknown", colnames(single_integrals))]

#rownames(single_integrals) = paste0(single_integrals$species,single_integrals$experiment_id)
rownames(single_integrals) = single_integrals$species
#single_integrals = single_integrals[single_integrals$species %in% single_species,]
# Remove everything but metabolites
mets_no_unknown = mets[-grep("unknown",mets)]
single_integrals = single_integrals[,mets_no_unknown]
# only keep potentially cross-fed metabolites
single_integrals[single_integrals > 10] = 10
single_integrals[single_integrals < -10] = -10
single_integrals[single_integrals < 1 & single_integrals > -1] = 0
single_integrals = single_integrals[,sapply(single_integrals,max) > 2 & sapply(single_integrals,min) < -2]

pheatmap(single_integrals,filename='~/last_metabolomics_heatmap.jpg',width=9,height=4,colorRampPalette(c("blue", "white", "red"))(100))



# For each species, calculate a yield per genome






plot(single_integrals[,mets[13]],log_all_data[log_all_data$species %in% single_species,]$normalized_OD)
mets[1:13]

# generate linear model for each species, predicting metabolite abundance using species OD
#log-transformed z-scores
#filter by samples with OD > 0.01
log_od_filtered = log_all_data[log_all_data$normalized_OD > 0.01,]
# replace 2-Oxoisocaproate with 2_Oxoisocaproate to avoid formula issue in linear model
colnames(log_od_filtered)[colnames(log_od_filtered) == '2-Oxoisocaproate'] = 'Oxoisocaproate'
mets[mets=='2-Oxoisocaproate'] = 'Oxoisocaproate'
integrals_356 = log_od_filtered[log_od_filtered$species %in% c('356'),]
integrals_360 = log_od_filtered[log_od_filtered$species %in% c('360'),]
integrals_361 = log_od_filtered[log_od_filtered$species %in% c('361'),]
integrals_492 = log_od_filtered[log_od_filtered$species %in% c('492'),]
integrals_500 = log_od_filtered[log_od_filtered$species %in% c('500'),]
integrals_502 = log_od_filtered[log_od_filtered$species %in% c('502'),]
integrals_519 = log_od_filtered[log_od_filtered$species %in% c('519'),]

# Generate a filtered dataframe for each species and calculate a yield
# using log-transformed (non-z) metabolite integral/DNA quantification value


# Generate new predicted log(metabolite integral values) from co-cultures using
# an additive assumption

# Z-transform the co-culture predictions and compare to observed

for (metabolite in mets) {
  print(metabolite)
  fmla <- paste0("normalized_OD ~ ",metabolite)
  linmod = lm(fmla, data = integrals_356)
}
fmla <- as.formula(paste("normalized_OD ~ ", paste(mets, collapse= " + ")))
fmla <- ("normalized_OD ~ unknown1")
form = lm(fmla, data = integrals_356)

integrals_356


