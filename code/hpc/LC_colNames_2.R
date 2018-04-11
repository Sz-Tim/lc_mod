# After running LC_colNames.sh to fill full_colNames/
# This identifies which columns to exclude
# There is a way to do this in the same bash script to avoid 
# the back and forth...

col.dir <- "full_colNames/"
exc.dir <- "full_colExc/"
inc.dir <- "full_colInc/"

m.i <- list.files(col.dir, pattern="^full_Clim_")

# identify columns to include in processing
for(i in 1:length(m.i)) {
  v.nm <- scan(paste0(col.dir, m.i[i]), sep=",", what="character")
  col.inc <- grep(paste(c("nu", "Y2", "Y2_p_", "Y2_"), collapse="|"), v.nm, invert=T)
  cat(col.inc, paste0(inc.dir, m.i[i]), sep=",")
}

# identify continuous sequences of indices for bash::cut in LC_colNames_3.sh
col.inc <- scan("full_colInc/full.csv", sep=",")
col.inc <- col.inc[-which(is.na(col.inc))]
lag.i <- col.inc[-1] - col.inc[-length(col.inc)]
which(lag.i > 1) # 1-67, 68-84, 85-1444032, 1444033-1444048 --- col.inc[]
col.vals <- list(c(col.inc[1], col.inc[67]),
                 c(col.inc[68], col.inc[84]),
                 c(col.inc[85], col.inc[1444032]),
                 c(col.inc[1444033], col.inc[length(col.inc)]))