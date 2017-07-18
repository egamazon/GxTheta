SimulateDifferentialTagging <-
function(save.path, n.iter, genetic.effect.size, n.idvs) {
    cols <- c('ait.chi2', 'aitl.chi2', 'r.eur', 'r.afr')
    results <- matrix(NA, nrow=n.iter, ncol=length(cols))
    colnames(results) <- cols
    for (n in 1:n.iter) {
        # define haplotypes and frequencies of each haplotype
        eur.h1.frq <- 0.4
        eur.h2.frq <- 0.5
        eur.h3.frq <- 0.1

        eur.h1 <- c(0, 0)
        eur.h2 <- c(0, 1)
        eur.h3 <- c(1, 1)
        eur.haps <- t(data.frame(eur.h1, eur.h2, eur.h3))

        afr.h1.frq <- 0.01
        afr.h2.frq <- 0.4
        afr.h3.frq <- 0.01

        afr.h1 <- c(1, 1)
        afr.h2 <- c(1, 0)
        afr.h3 <- c(0, 0)
        afr.haps <- t(data.frame(afr.h1, afr.h2, afr.h3))

        # draw global ancestry
        global.euro.anc <- rnorm(n.idvs, 0.7, 0.2)
        global.euro.anc[ global.euro.anc > 1 ] <- 1
        global.euro.anc[ global.euro.anc < 0 ] <- 0

        # draw haplotypes for each individual
        local.euro.anc.h1 <- rbinom(n.idvs, 1, global.euro.anc)
        local.euro.anc.h2 <- rbinom(n.idvs, 1, global.euro.anc)

        local.euro.anc <- local.euro.anc.h1 + local.euro.anc.h2

        n.eur.haps <- sum(local.euro.anc)
        n.afr.haps <- (2*n.idvs)-n.eur.haps

        eur.haps.check <- matrix(NA, nrow=n.eur.haps, ncol=2) 
        afr.haps.check <- matrix(NA, nrow=n.afr.haps, ncol=2) 

        genos <- matrix(NA, nrow=n.idvs, ncol=2)
        genos.idx <- 1
        eur.idx <- 1
        afr.idx <- 1
        for (i in local.euro.anc) {
            if (i == 2) {
                # draw from european haplotypes
                h1 <- which(rmultinom(1, 1, c(eur.h1.frq, eur.h2.frq, eur.h3.frq)) == 1)
                h2 <- which(rmultinom(1, 1, c(eur.h1.frq, eur.h2.frq, eur.h3.frq)) == 1)
                hap_genos <- eur.haps[h1, ] + eur.haps[h2, ]
                genos[genos.idx, ] <- hap_genos
                eur.haps.check[eur.idx, ] <- eur.haps[h1, ]
                eur.idx <- eur.idx + 1
                eur.haps.check[eur.idx, ] <- eur.haps[h2, ]
                eur.idx <- eur.idx + 1
            }
            else if (i==1) {
                # draw from african and european
                h1 <- which(rmultinom(1, 1, c(eur.h1.frq, eur.h2.frq, eur.h3.frq)) == 1)
                h2 <- which(rmultinom(1, 1, c(afr.h1.frq, afr.h2.frq, afr.h3.frq)) == 1)
                hap_genos <- eur.haps[h1, ] + afr.haps[h2, ]
                genos[genos.idx, ] <- hap_genos
                eur.haps.check[eur.idx, ] <- eur.haps[h1, ]
                eur.idx <- eur.idx + 1
                afr.haps.check[afr.idx, ] <- afr.haps[h2, ]
                afr.idx <- afr.idx + 1
            }
            else {
                # draw from an african haplotypes
                h1 <- which(rmultinom(1, 1, c(afr.h1.frq, afr.h2.frq, afr.h3.frq)) == 1)
                h2 <- which(rmultinom(1, 1, c(afr.h1.frq, afr.h2.frq, afr.h3.frq)) == 1)
                hap_genos <- afr.haps[h1, ] + afr.haps[h2, ]
                genos[genos.idx, ] <- hap_genos
                afr.haps.check[afr.idx, ] <- afr.haps[h1, ]
                afr.idx <- afr.idx + 1
                afr.haps.check[afr.idx, ] <- afr.haps[h2, ]
                afr.idx <- afr.idx + 1
            }

            genos.idx <- genos.idx + 1
        }

        # check what kind of correlations we are getting
        r.eur <- cor(eur.haps.check[,1], eur.haps.check[,2])
        r.afr <- cor(afr.haps.check[,1], afr.haps.check[,2])

        maf.eur.g1 <- sum(eur.haps.check[,1])/(nrow(eur.haps.check))
        if (maf.eur.g1 > 0.5) {
            maf.eur.g1 <- 1-maf.eur.g1
        }
        maf.eur.g2 <- sum(eur.haps.check[,2])/(nrow(eur.haps.check))
        if (maf.eur.g2 > 0.5) {
            maf.eur.g2 <- 1-maf.eur.g2
        }
        maf.afr.g1 <- sum(afr.haps.check[,1])/(nrow(afr.haps.check))
        if (maf.afr.g1 > 0.5) {
            maf.afr.g1 <- 1-maf.afr.g1
        }
        maf.afr.g2 <- sum(afr.haps.check[,2])/(nrow(afr.haps.check))
        if (maf.afr.g2 > 0.5) {
            maf.afr.g2 <- 1-maf.afr.g2
        }

        # generate some phenotypes based on the 2nd genotype of each individual
        phenos <- rnorm(n.idvs, mean=genetic.effect.size*genos[,2], sd=1)

        # create data frame
        idv.data <- data.frame(phenos, genos, global.euro.anc, local.euro.anc)

        # do the testing
        test1 <- summary(lm(formula="phenos ~ X1 + global.euro.anc + X1:global.euro.anc", data=idv.data))
        test2 <- summary(lm(formula="phenos ~ X1 + global.euro.anc + X1:global.euro.anc + local.euro.anc", data=idv.data))

        p.g1 <- (maf.eur.g1 + maf.afr.g1)/2
        fst.g1 <- (maf.eur.g1 - maf.afr.g1)^2 / (2*p.g1*(1-p.g1))
        p.g2 <- (maf.eur.g2 + maf.afr.g2)/2
        fst.g2 <- (maf.eur.g2 - maf.afr.g2)^2 / (2*p.g2*(1-p.g2))

        test.results <- c((test1$coefficients[4,3])^2, (test2$coefficients[5,3])^2, r.eur, r.afr)
        results[n, ] <- test.results
    }

    write.table(results, save.path, row.names=F, col.names=T, quote=F)
}
