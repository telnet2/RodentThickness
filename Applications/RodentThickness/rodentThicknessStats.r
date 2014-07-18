########################################################
# Script to compute population thickness difference
#
# Rscript rodentThicknessStats.r work-dir control-data treatment-data
options(warn=-1)
cmdArgs = commandArgs(TRUE);
ctlFile = cmdArgs[1]; 
expFile = cmdArgs[2];
outfile = cmdArgs[3];
firstgroup = cmdArgs[4];
secondgroup = cmdArgs[5];
cl <- data.matrix(read.csv(file=ctlFile, sep="\t", header=FALSE));
fl <- data.matrix(read.csv(file=expFile, sep="\t", header=FALSE));

# remove the trailing delimiter
#cl = cl[,1:(ncol(cl)-1)];
#fl = fl[,1:(ncol(fl)-1)];

cl = cl[,1:(ncol(cl))];
fl = fl[,1:(ncol(fl))];

alldata = cbind(cl,fl);
clCols = 1:ncol(cl);
flCols = (ncol(cl)+1):(ncol(fl) + ncol(cl));

#########################################################
# compute each group's row-wise (point-wise) mean
print("Start t-test")
alltest = apply(alldata,1,function(x) t.test(x[clCols],x[flCols]));

pvalues = sapply(alltest, function(x) x$p.value);
tscores = sapply(alltest, function(x) x$statistic);
adjpvalues = p.adjust(pvalues, method='fdr');

clMeans = sapply(alltest, function(x) x$estimate[1]);
flMeans = sapply(alltest, function(x) x$estimate[2]);

print("Start wilcox-test")
ranktest = apply(alldata,1,function(x) tryCatch(wilcox.test(x[clCols],x[flCols]), error = 1));
print("done..")

wpvalues = sapply(ranktest, function(x) x$p.value);
wscores = sapply(ranktest, function(x) x$statistic);
print("FDR ...")
wadjpvalues = p.adjust(wpvalues, method='fdr');
print("done ...")


flStd = apply(fl,1,sd);
clStd = apply(cl,1,sd);

meanDiff = clMeans - flMeans;
stdDiff  = clStd - flStd;

print(pvalues)
filter_significant <- function(x) {
  if (is.nan(x)) {
    return(NaN);
  } else if (x < 0.05) {
    return(1);
  } else {
    return(0);
  }
}
clRatio = 100*meanDiff / clMeans * sapply(pvalues, filter_significant);
flRatio = 100*meanDiff / flMeans * sapply(pvalues, filter_significant);

meanDiffMask = meanDiff * sapply(pvalues, filter_significant);
meanDiffIndex = sapply(meanDiffMask, function(x) if (is.nan(x)) { return(NaN); } else if (x >= 0) { floor(x / 100) * 100 } else { ceiling(x / 100) * 100 });


name1 = paste(c(firstgroup, "Mean"), collapse = "_") ;
name2 = paste(c(secondgroup, "Mean"), collapse = "_") ;
name3 = paste(c(name1, name2), collapse = "-") ;
name4 = paste(c(firstgroup, "Sd"), collapse = "_") ;
name5 = paste(c(secondgroup, "Sd"), collapse = "_") ;
name6 = paste(c(firstgroup, "Ratio"), collapse = "_");
name7 = paste(c(secondgroup, "Ratio"), collapse = "_");
name8 = paste(c(firstgroup, secondgroup, "Index"), collapse = "_");

# manual t-test computation
#ttest = (clMean-flMean)/(clStd^2/ncol(cl)+flStd^2/ncol(fl))^0.5;
#pvalues = pt(ttest, ncol(cl)+ncol(fl)-2);

outputTable = cbind(clMeans,flMeans,meanDiff,clStd,flStd,tscores,pvalues,adjpvalues,wscores,wpvalues,wadjpvalues, clRatio, flRatio, meanDiffIndex);
colnames(outputTable) = c(name1, name2, name3, name4, name5, "t.scores", "t.pvalue", "t.adjpvalue", "w.scores", "w.pvalue", "w.adjpvalue", name6, name7, name8);
write.csv(outputTable, row.names=FALSE, file=outfile, na="NaN");
