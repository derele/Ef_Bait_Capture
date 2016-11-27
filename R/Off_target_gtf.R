

baits <- read.table("/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_feature_counted.gtf")


test.case <- baits[baits$V1%in%"EfaB_7715",]

add.off.target <- function(x) {
    ## order by the start of the first feature
    ordered <- x[order(x$V4), ]
    ordered
}



by(test.case, as.character(test.case2$V1), function (x) add.off.target(x))

test.case2 <- baits[baits$V1%in%c("EfaB_7715","EfaB_6072"),]

by(test.case2, as.character(test.case2$V1), function (x) add.off.target(x))
