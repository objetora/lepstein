Sys.setlocale("LC_MESSAGES", "C")

# define the functions
wilcox.p = function (x) wilcox.test(x)[3] # extract wilcoxon p-value
nonsig = function(x) x >= 0.05 # TRUE or FALSE if significant

# load the tables
# pcri has the initial extraction values
pcri = read.table("inicial.csv", header=TRUE, sep="\t")
# varin has the pre values calculated from the initial extraction
# values and the mixture percentages, use it to compute
# multipliers to convert extractions into mixture percents
varin = read.table("varin.csv", header=TRUE, sep="\t")
# pcr has the normalized cpt values for each method, mixture, rep, and isolate
pcr = read.table("pcr-norm.csv", header=TRUE, sep="\t")

# differences for observed data
pcr$diff = pcr$post - pcr$pre

# set up column
pcr$methisol = paste("y", pcr$meth, pcr$isol, sep="_")
obsmedians = data.frame(unlist(tapply(pcr$diff, pcr$methisol, median))) / 10
obsmedians$methisol = rownames(obsmedians)
obsmedians$meth = substr(obsmedians$methisol, 3, 6)
obsmedians$isol = substr(obsmedians$methisol, 8, 100)
colnames(obsmedians) = c("median", "methisol", "meth", "isol")
obsq13s = data.frame(unlist(tapply(pcr$diff, pcr$methisol, quantile, probs = .25))) / 10
obsq13s$methisol = rownames(obsq13s)
obsq13s$meth = substr(obsq13s$methisol, 3, 6)
obsq13s$isol = substr(obsq13s$methisol, 8, 100)
colnames(obsq13s) = c("q1", "methisol", "meth", "isol")
obsq3s = data.frame(unlist(tapply(pcr$diff, pcr$methisol, quantile, probs = .75))) / 10
colnames(obsq3s) = "q3"
obsq13s$q3 = obsq3s


# compute wilcoxon, adjusted for BH for all isolates, all methods
obswilcoxvals = t(data.frame(tapply(pcr$diff, pcr$methisol, wilcox.p)))
obswilcoxs = data.frame(cbind(rownames(obswilcoxvals), obswilcoxvals))
colnames(obswilcoxs) = c("methisol", "wilcox")
obswilcoxs$meth = substr(obswilcoxs$methisol, 1, 6)
obswilcoxs$isol = substr(obswilcoxs$methisol, 7, 100)
obswilcoxs$wilcox = as.numeric(obswilcoxs$wilcox)
dim(obswilcoxs$wilcox) = nrow(obswilcoxs)
obswilcoxsbh = data.frame(unlist(p.adjust(obswilcoxs$wilcox, method="BH")))
obswilcoxsbh$methisol = obswilcoxs$methisol
obswilcoxsbh$isol = obswilcoxs$isol
obswilcoxsbh$meth = obswilcoxs$meth
colnames(obswilcoxsbh) = c("wilcoxon_BH", "methisol", "isol", "meth")
obswilcoxsbh$methisol = obswilcoxs$methisol
obswilcoxsbh$meth = obswilcoxs$meth
obswilcoxsbh$isol = obswilcoxs$isol
obswilcoxs = obswilcoxs[order(obswilcoxs$isol, obswilcoxs$meth), ]
obswilcoxsbh = obswilcoxsbh[order(obswilcoxsbh$isol, obswilcoxsbh$meth), ]

# prepare to compute multipliers to convert extractions
# into mixture percents
pcr$isolmethtmt = paste(pcr$isol, pcr$meth, pcr$tmt, sep="_")
pcri$isolrep = paste(pcri$isol, pcri$rep, sep="_")
varin$isolrep = paste(varin$isol, varin$rep, sep="_")
pcrv = merge(pcri, varin, by=c("isol", "rep", "isolrep"))
pcrv$rat = pcrv$pre / pcrv$extrac
pcrv$isolmix = paste(pcrv$isol, pcrv$mix, sep="_")

# prems are the multipliers
premns = data.frame(tapply(pcrv$pre, pcrv$isolmix, mean))
isolmix = rownames(premns)
premns$isolmix = isolmix
colnames(premns) = c("premn", "isolmix")
pcrv = merge(pcrv, premns, by="isolmix")
pcrv$emrat = pcrv$rat / pcrv$premn
emrat = data.frame(pcrv[pcrv$mix == 1, c("isol", "extrac", "emrat")])

# check normality of emrats
sw = tapply(emrat$extrac, emrat$isol, shapiro.test)

# create initialization data frame for the normal distributions
# to sample from
mn = data.frame(tapply(emrat$extrac, emrat$isol, mean))
mn$isol = rownames(mn)
colnames(mn) = cbind("mn", "isol")
sd = as.data.frame(tapply(emrat$extrac, emrat$isol, sd))
sd$isol = rownames(sd)
colnames(sd) = cbind("sd", "isol")
inic = as.data.frame(merge(emrat, mn, by="isol"))
inic = as.data.frame(merge(inic, sd, by="isol"))
inic = as.data.frame(inic[!duplicated(inic[ ,
													cbind("isol", "mn")]), colnames(inic) != "extrac"])
for (s in c(1:10000))
{
	# create random extractions
	inic$ext = mapply(rnorm, n = 1, mean = inic$mn, sd = inic$sd, SIMPLIFY=TRUE)
	pcrrand = merge(pcr, inic[ , cbind("isol", "ext", "emrat")], by="isol")

	# convert to pre values
	pcrrand$pre = pcrrand$ext * pcrrand$emrat * pcrrand$pre

	#
	pcrrand$methtmt = paste("y", pcrrand$meth, pcrrand$tmt, sep="_")
	sumpre = data.frame(tapply(pcrrand$pre, pcrrand$methtmt, sum))
	sumpre$methtmt = rownames(sumpre)
	colnames(sumpre) = cbind("sumpre", "methtmt")
	pcrrand = data.frame(merge(pcrrand, sumpre, by="methtmt"))
	# scale to 1000 for cpt
	pcrrand$pre = pcrrand$pre * 1000 / pcrrand$sumpre

	# compute wilcoxon, adjusted by BH, with random pre
	pcrrand$diff = pcrrand$post - pcrrand$pre
	wilcoxvals = t(data.frame(tapply(pcrrand$diff, pcrrand$methisol, wilcox.p)))
	wilcoxs = data.frame(cbind(rownames(wilcoxvals), wilcoxvals))
	colnames(wilcoxs) = c("methisol", "wilcox")
	wilcoxs$meth = substr(wilcoxs$methisol, 1, 6)
	wilcoxs$isol = substr(wilcoxs$methisol, 7, 100)
	wilcoxs$wilcox = as.numeric(wilcoxs$wilcox)
	wilcoxsbh = data.frame(unlist(p.adjust(wilcoxs$wilcox, method="BH")), row.names=wilcoxs$methisol)
	colnames(wilcoxsbh) = "wilcoxon_BH"
	wilcoxsbh$methisol = wilcoxs$methisol
	wilcoxsbh$meth = wilcoxs$meth
	wilcoxsbh$isol = wilcoxs$isol

	# add to the list of means of wilcoxs
	if (exists("wilcoxrand"))
	{
		wilcoxrand1 = data.frame(t(wilcoxsbh$wilcox), row.names=s)
		colnames(wilcoxrand1) = wilcoxsbh$methisol
		wilcoxrand = rbind(wilcoxrand, wilcoxrand1)
	}
	else
	{
		wilcoxrand = data.frame(t(wilcoxsbh$wilcox), row.names=1)
		colnames(wilcoxrand) = wilcoxsbh$methisol
	}
}
methisol = colnames(wilcoxrand)
wilcoxrandnonsig = nonsig(wilcoxrand)
wilcoxrandnonsigfract = data.frame(unlist(apply(wilcoxrandnonsig, 2, mean)))
colnames(wilcoxrandnonsigfract) = "fractnonsig"
wilcoxrandnonsigfract$methisol = methisol
wilcoxrandnonsigfract$isol = wilcoxsbh$isol
wilcoxrandnonsigfract$meth = wilcoxsbh$meth
