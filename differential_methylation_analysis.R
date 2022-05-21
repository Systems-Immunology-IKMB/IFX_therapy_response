library(RnBeads)
library(LOLA)
library(GO.db)

#Code for differential methylation analysis
#Parameters can be modified for different timepoints and disease status groups

#Specify project directory
data.dir <- "Epicarray/"

#Specify directory containing iDat files
idat.dir <- file.path(data.dir, "iDat_files/")

#Specify analysis directory 
analysis.dir <- file.path(data.dir, "Analysis/IFX/comparison_0h_2w/R/")

#Specify report directories for each comparison
report.dir <- file.path(analysis.dir, "reports_differential/")

#Specify annotation files for each comparison
sample.annotation <- file.path(data.dir, "EpicArray_IFX_R_0h_2w.csv")
data.type<-"idat.dir"

#Set rnbeads options

help("rnb.options")
rnb.options(analysis.name = "example")
rnb.options(logging = TRUE)
rnb.options(assembly = "hg19")
rnb.options(identifiers.column = "Sample_ID")
rnb.options2xml(pretty=TRUE)
rnb.options(min.group.size = 1)

###### modules
rnb.options(import = TRUE)
rnb.options(preprocessing = TRUE)
rnb.options(qc = TRUE)
rnb.options(inference = FALSE)
rnb.options(exploratory = FALSE)
rnb.options(differential = TRUE)

##############some selected options (see more at help("rnb.options"))

rnb.options(import.default.data.type = "infinium.idat.dir") # NGS: "bs.bed.dir"
rnb.options(import.table.separator = ",") # for TABdel: "\t"
rnb.options(import.gender.prediction = FALSE) # only for array datasets

#normalization
rnb.options(normalization = NULL)
rnb.options(normalization.method = "illumina")
rnb.options(normalization.background.method = "methylumi.noob")
rnb.options(normalization.plot.shifts = TRUE)

#qc
rnb.options(qc.boxplots = TRUE)
rnb.options(qc.barplots = TRUE)
rnb.options(qc.negative.boxplot = TRUE)
rnb.options(qc.snp.distances = TRUE)
rnb.options(qc.snp.boxplot = TRUE)
rnb.options(qc.snp.barplot = TRUE)
rnb.options(qc.sample.batch.size = 50)
rnb.options(qc.coverage.plots = FALSE)
rnb.options(qc.coverage.threshold.plot = 1:10)
rnb.options(qc.coverage.histograms = FALSE)
rnb.options(qc.coverage.violins = FALSE)

#filtering
rnb.options(filtering.whitelist = NULL)
rnb.options(filtering.blacklist = NULL)
rnb.options(filtering.snp = "3")
rnb.options(filtering.cross.reactive = FALSE)
rnb.options(filtering.greedycut = TRUE)
rnb.options(filtering.greedycut.pvalue.threshold = 0.05)
rnb.options(filtering.greedycut.rc.ties = "row")
rnb.options(filtering.sex.chromosomes.removal = TRUE)
rnb.options(filtering.missing.value.quantile = 0.8)
rnb.options(filtering.coverage.threshold = 5)
rnb.options(filtering.low.coverage.masking = FALSE)
rnb.options(filtering.high.coverage.outliers = FALSE)
rnb.options(filtering.deviation.threshold = 0)

#regions
rnb.load.annotation.from.db("ensembleRegBuildBPall")
rnb.load.annotation.from.db("ensembleRegBuildBPproximal")
rnb.load.annotation.from.db("ensembleRegBuildBPdistal")
rnb.options(region.types=c("promoters","genes","cpgislands","tiling","ensembleRegBuildBPall","ensembleRegBuildBPproximal","ensembleRegBuildBPdistal")) 

#differential analysis
rnb.options(differential.site.test.method = "ttest") 
rnb.options(differential.comparison.columns = c("Timepoint")) #Set column used for differential comparision

rnb.options(covariate.adjustment.columns = NULL)
rnb.options(columns.pairing = c("Timepoint"="Patient")) #For paired analyis
rnb.options(differential.adjustment.sva = TRUE)
rnb.options(differential.adjustment.celltype = FALSE)
rnb.options(differential.enrichment = FALSE)
rnb.options(differential.report.sites = TRUE)

rnb.options(export.to.bed = TRUE)
rnb.options(export.to.trackhub = c("bigBed","bigWig"))
rnb.options(export.to.csv = TRUE) #default = FALSE #creates big tables
rnb.options(export.to.ewasher = FALSE)
rnb.options(export.types = "sites")

rnb.options(disk.dump.big.matrices = TRUE)
rnb.options(enforce.memory.management = FALSE)
rnb.options(enforce.destroy.disk.dumps = FALSE)

##############################################################################

##### Run Analysis
rnb.run.analysis(dir.reports=report.dir,
                 sample.sheet=sample.annotation,
                 data.dir=idat.dir,
                 data.type=data.type)


##############################################################################


