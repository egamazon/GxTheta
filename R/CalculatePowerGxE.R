CalculatePowerGxE <-
function (save.path, pval.cutoff, E.effect.size, genetic.effect.size, global.anc.effect, interaction.effect.size, n.idvs, sigma.max, sigma.step, n.iter) {
    rows <- n.iter * (sigma.max / sigma.step) + n.iter
    cols <- c('direct.test1.chi2', 'direct.test2.chi2', 'test1.chi2', 'test2.chi2', 'test1.pval', 'test2.pval', 'r.etheta', 'sigma', 'p.a', 'p.e', 'p.ae')
    results <- matrix(NA, nrow=rows, ncol=length(cols))
    results.idx <- 1

    for (n in 1:n.iter) {
        print(n)
        # set global ancestry
        global.euro.anc <- rnorm(n.idvs, 0.7, 0.2)
        global.euro.anc[ global.euro.anc > 1 ] <- 1
        global.euro.anc[ global.euro.anc < 0 ] <- 0

        # set local ancestry
        local.euro.anc.h1 <- matrix(rbinom(n.idvs, 1, global.euro.anc), nrow=n.idvs, ncol=1)
        local.euro.anc.h2 <- matrix(rbinom(n.idvs, 1, global.euro.anc), nrow=n.idvs, ncol=1)
        local.euro.anc <- local.euro.anc.h1 + local.euro.anc.h2

        # set allele frequencies for ancestries
        mafs1.euro = mafs1.afr = 0
        p.ae <- 0
        p.e <- 0
        p.a <- 0
        Fst <- 0.16

        while (p.e < 0.05 || p.e > 0.95) {
            p.e <- rbeta(1, 0.2*(1-Fst)/Fst,(1-0.2)*(1-Fst)/Fst)
        }

        while (p.a < 0.05 || p.a > 0.95) {
            p.a <- rbeta(1, 0.2*(1-Fst)/Fst,(1-0.2)*(1-Fst)/Fst)
        }

        p.ae <- (p.a + p.e) / 2

        # set genotypes
        mafs1.euro[1] <- p.e
        mafs1.afr[1] <- p.a

        # set genotypes for the first SNP

        h1 <- matrix(0, nrow=n.idvs, ncol=1);
        h2 <- matrix(0, nrow=n.idvs, ncol=1);
        for(i in 1:n.idvs) {
            #haplotype maternal
            h1[i,local.euro.anc.h1[i, 1]==0] <- rbinom(sum(local.euro.anc.h1[i, 1]==0), 1, mafs1.afr[local.euro.anc.h1[i, 1]==0]) 
            h1[i,local.euro.anc.h1[i, 1]==1] <- rbinom(sum(local.euro.anc.h1[i, 1]==1), 1, mafs1.euro[local.euro.anc.h1[i, 1]==1])

            #haplotype paternal
            h2[i,local.euro.anc.h2[i, 1]==0] <- rbinom(sum(local.euro.anc.h2[i, 1]==0), 1, mafs1.afr[local.euro.anc.h2[i, 1]==0])
            h2[i,local.euro.anc.h2[i, 1]==1] <- rbinom(sum(local.euro.anc.h2[i, 1]==1), 1, mafs1.euro[local.euro.anc.h2[i, 1]==1])
        }

        genos <- data.frame(h1+h2)

        # generate enviromental exposures
        for (s in seq(0, sigma.max, sigma.step)) {
            E <- scale(rnorm(n.idvs, mean=global.euro.anc, sd=s))
            r.etheta <- cor(global.euro.anc, E)

            phenos.mean <- (genos[,1]*genetic.effect.size) + (global.euro.anc * global.anc.effect) + (E * E.effect.size) + (genos[,1] * E * interaction.effect.size)
            phenos <- rnorm(n.idvs, mean=phenos.mean, sd=1)

            g <- genos[, 1]
            idv.data <- data.frame(phenos, g, global.euro.anc, E)
            
            test <- summary(lm(formula="phenos ~ g + global.euro.anc + g:global.euro.anc", data=idv.data))
            direct.test <- summary(lm(formula="phenos ~ g + global.euro.anc + g:E", data=idv.data))

            test1 <- test$coefficients[4,3]^2 # this is 1 dof
            test1.pval <- 1-pchisq(test1, 1)

            test2 <- test$coefficients[4,3]^2 + test$coefficients[2,3]^2 # this is 2 dof
            test2.pval <- 1-pchisq(test2, 2)

            direct.test1 <- direct.test$coefficients[4,3]^2 # this is 1 dof
            direct.test2 <- direct.test$coefficients[4,3]^2 + direct.test$coefficients[2,3]^2 # this is 2 dof

            test.results <- c(direct.test1, direct.test2, test1, test2, test1.pval, test2.pval, r.etheta, s, p.a, p.e, p.ae)

            results[results.idx, ] <- test.results
            results.idx <- results.idx + 1
        }
    }

    colnames(results) <- cols
    results <- as.data.frame(results)

    power.results <- matrix(NA, nrow=sigma.max/sigma.step, ncol=3)
    chi2.cutoff <- qchisq(-1*(pval.cutoff-1),1)
    idx <- 1
    # # plot simulmation for 1-dof test
    for (i in seq(0, 1, 0.05)) {
        tmp <- results[ results$r.etheta >= i -  0.05 & results$r.etheta < i, ]
        tmp <- tmp[, c('test1.chi2', 'direct.test1.chi2', 'r.etheta')]
        power.test1 <- nrow(tmp[ tmp$test1.chi2 > chi2.cutoff, ]) / nrow(tmp)
        power.test2 <- nrow(tmp[ tmp$direct.test1.chi2 > chi2.cutoff, ]) / nrow(tmp)
        power.results[idx, ] <- c(power.test1, power.test2, i)
        idx <- idx + 1
    }

    power.results <- as.data.frame(power.results)
    colnames(power.results) <- c('ait.test1.power', 'direct.test1.power', 'r.etheta')
    power.results <- power.results[ !is.na(power.results$r.etheta), ]
    write.table(power.results, save.path, row.names=F, col.names=T, quote=F)
}
