# speed4
# Type of the joint rhythmicity -------------------------------------------
rm(list = ls())
hdir = "/home/xix66/circadian/ThePipeline/CircadianPipeline_private/AW_fisher" #speed4
# dir.create(hdir, recursive = TRUE)
out.dir = paste0(hdir, "/output_TOJRsim")
dir.create(out.dir, recursive = TRUE)
ncores = 30
# install.packages("devtools")    # if it is not already installed
# devtools::install_github("cran/npsm")   # Package archived by CRAN
# devtools::install_github("cran/DODR")   # Package archived by CRAN
# BiocManager::install("rain")
# BiocManager::install("DESeq2")
# devtools::install_github("bharathananth/compareRhythms")
# source(paste0(hdir, "/code/CP_DiffRhythmicity.R"))
# source(paste0(hdir, "/code/New packages.R"))
# source(paste0(hdir, "/code/two_cosinor_OLS_overall.R"))
# devtools::install_github("diffCircadian/diffCircadian")
# devtools::install_github("Caleb-Huo/differentialR2")
# devtools::install_github("Caleb-Huo/AWFisher")
# library(diffCircadian)
# library(differentialR2)
# install.packages("selectiveInference")
library(compareRhythms)
library(tidyverse)
library(parallel)
source(paste0(hdir, "/code/compareRhythm_dodr_replace.R"))

# Settings for types of the joint rhythmicity -----------------------------
# 1. arhythmic: A1 = A2 = 0, sigma1 = sigma2 = 1, M1 = M2 = 5
# 2. rhythmic in group I: A1 = 3; A2 = 0; sigma1 = sigma2 = 1; phase1 = 0; M1 = M2 = 5
# 3. rhythmic in group II: A1 = 0; A2 = 3; sigma1 = sigma2 = 1; phase2 = 0; M1 = M2 = 5
# 4. rhythmic in both with no change of rhythmicity: A1 = 3; A2 = 3; sigma1 = sigma2 = 1; phase1 = phase2 = 0; M1 = M2 = 5
# 5. rhythmic in both with change of parameters: A1 = 3; A2 = 3; sigma1 = sigma2 = 1; phase1 = 0; phase2 = pi/2, ...; M1 = M2 = 5
  # 5. rhythmic in both with change of parameters: A1 = 3; A2 = 3; sigma1 = sigma2 = 1; phase1 = 0; phase2 = pi/2, ...; M1 = M2 = 5
# 6. rhythmic in both with change of fitness: A1 = 3; A2 = 3; sigma1 = 1.5; sigma2 = 1; phase1 = 0; phase2 = pi/2; M1 = M2 = 5

categories.all = c("arrhy",
                   "rhy1.0.5", "rhy1.0.6", "rhy1.0.7", "rhy1.0.8", "rhy1.0.9", "rhy1.1",
                   "same.0.5", "same.0.6", "same.0.7", "same.0.8", "same.0.9", "same.1",
                   "diffPar.4", "diffPar.6", "diffPar.8", "diffPar.10",
                   "diffR2.0.67", "diffR2.0.83", "diffR2.1.25", "diffR2.1.67", "diffR2.2.50")

settings = lapply(categories.all, function(a.cat){
  if(a.cat == "arrhy"){
    a.row = data.frame(A1 = 0, A2  = 0, sigma1 = 1, sigma2 = 1, phase1 = 0, phase2 = 0, M1 = 5, M2 = 5, setting = a.cat)
  }else if(grepl("rhy1", a.cat)){
    a.row = data.frame(A1 = as.numeric(gsub("rhy1.([0-9\\.]+)$", "\\1", a.cat)),
                       A2  = 0, sigma1 = 1, sigma2 = 1, phase1 = 0, phase2 = 0, M1 = 5, M2 = 5, setting = a.cat)
  }else if(grepl("same", a.cat)){
    a.A = as.numeric(gsub("same.([0-9\\.]+)$", "\\1", a.cat))
    a.row = data.frame(A1 = a.A, A2  = a.A, sigma1 = 1, sigma2 = 1, phase1 = 0, phase2 = 0, M1 = 5, M2 = 5, setting = a.cat)
  }else if(grepl("diffPar", a.cat)){
    a.phase = as.numeric(gsub("diffPar.([0-9\\.]+)$", "\\1", a.cat))
    a.row = data.frame(A1 = 1, A2  = 1, sigma1 = 1, sigma2 = 1, phase1 = pi/a.phase, phase2 = 0, M1 = 5, M2 = 5, setting = a.cat)
  }else if(grepl("diffR2", a.cat)){
    a.sigma = as.numeric(gsub("diffR2.([0-9\\.]+)$", "\\1", a.cat))
    a.row = data.frame(A1 = 1, A2  = 1, sigma1 = a.sigma, sigma2 = 1, phase1 = 0, phase2 = 0, M1 = 5, M2 = 5, setting = a.cat)
  }
  return(a.row)
})
settings = do.call(rbind.data.frame, settings)

num_genes = 10000*10
nn = 30
tt = rep(c(0, 4, 8, 12, 16, 20), each = 5)

time.start = Sys.time()
res = mclapply(1:nrow(settings), function(a){
  set.seed(a)
  if(file.exists(paste0(out.dir, "/sim_cate", a, ".rds"))){
    NULL
  }else{
    # a.setting = lapply(1, function(a){
    a.category = settings$setting[a]
    a.A1 = settings$A1[a]
    a.A2 = settings$A2[a]
    a.sigma1 = settings$sigma1[a]
    a.sigma2 = settings$sigma2[a]
    a.phase1 = settings$phase1[a]
    a.phase2 = settings$phase2[a]
    a.M1 = runif(num_genes, min = 5, max = 10)
    a.M2 = runif(num_genes, min = 5, max = 10)

    a.dat1 = do.call(rbind, lapply(1:num_genes, function(i){
      a.A1*cos(2*pi/24*(tt-a.phase1))+a.M1[i]+rnorm(nn, 0, a.sigma1)
    }))
    rownames(a.dat1) = paste0(a.category,"_", seq_len(num_genes))
    a.dat2 = do.call(rbind, lapply(1:num_genes, function(i){
      a.A2*cos(2*pi/24*(tt-a.phase1))+a.M2[i]+rnorm(nn, 0, a.sigma2)
    }))
    rownames(a.dat2) = paste0(a.category,"_", seq_len(num_genes))

    a.dat = cbind(a.dat1, a.dat2)
    # tests -------------------------------------------------------------------
    #Sidak_BS, Sidak_FS, VDA, AWFisher
    x1 = list(data = data.frame(a.dat1),
              time = tt,
              gname = paste("gene", seq_len(nrow(a.dat1))))
    x2 = list(data = data.frame(a.dat2),
              time = tt,
              gname = paste("gene", seq_len(nrow(a.dat2))))

    aw.type = DiffCircaPipeline::DCP_Rhythmicity(x1, x2, method = "AWFisher", period = 24, alpha = 0.05, CI = FALSE, p.adjust.method = "BH")
    vda.type = DiffCircaPipeline::DCP_Rhythmicity(x1, x2, method = "VDA", period = 24, alpha = 0.05, CI = FALSE, p.adjust.method = "BH")
    Sidak_BS = DiffCircaPipeline::DCP_Rhythmicity(x1, x2, method = "Sidak_BS", period = 24, alpha = 0.05, CI = FALSE, p.adjust.method = "BH")
    Sidak_FS = DiffCircaPipeline::DCP_Rhythmicity(x1, x2, method = "Sidak_FS", period = 24, alpha = 0.05, CI = FALSE, p.adjust.method = "BH")

    a.meta = data.frame(time = c(tt, tt), group = rep(c("gI", "gII"), each = nn))
    results.1.bic <- compareRhythms(a.dat, exp_design = a.meta, amp_cutoff = 0,schwarz_wt_cutoff = 0,criterion = "bic",
                                    period = 24, method = "mod_sel",just_classify = FALSE) #changed the schwarz_wt_cutoff to 0. Otherwise not all genes are classified
    results.1.bic$category = as.character(results.1.bic$category)
    results.1.bic$category[results.1.bic$category=="loss"]="rhyI"
    results.1.bic$category[results.1.bic$category=="gain"]="rhyII"

    results.1.aic <- compareRhythms(a.dat, exp_design = a.meta, amp_cutoff = 0,schwarz_wt_cutoff = 0,criterion = "aic",
                                    period = 24, method = "mod_sel",just_classify = FALSE) #changed the schwarz_wt_cutoff to 0. Otherwise not all genes are classified
    results.1.aic$category = as.character(results.1.aic$category)
    results.1.aic$category[results.1.aic$category=="loss"]="rhyI"
    results.1.aic$category[results.1.aic$category=="gain"]="rhyII"

    results.2 <- compareRhythms2(a.dat, exp_design = a.meta,amp_cutoff = 0,
                                period = 24, method = "dodr",just_classify = FALSE)
    results.2$category = as.character(results.2$category)
    results.2$category[results.2$category=="loss"]="rhyI"
    results.2$category[results.2$category=="gain"]="rhyII"

    a.res = data.frame(id = rownames(a.dat1),
                       aw.type = aw.type$rhythm.joint$TOJR,
                       vda.type = vda.type$rhythm.joint$TOJR,
                       Sidak_BS = Sidak_BS$rhythm.joint$TOJR,
                       Sidak_FS = Sidak_FS$rhythm.joint$TOJR)
    colnames(results.1.bic)[2] = "mod_sel_bic"
    colnames(results.1.aic)[2] = "mod_sel_aic"
    colnames(results.2)[2] = "dodr"
    a.res = dplyr::left_join(a.res,  results.1.bic[, 1:2])
    a.res = dplyr::left_join(a.res,  results.1.aic[, 1:2])
    a.res = dplyr::left_join(a.res, results.2[, 1:2])
    a.res$truth = a.category
    a.res$dodr[is.na(a.res$dodr)] = "arrhy"
    saveRDS(a.res, paste0(out.dir, "/sim_cate", a, ".rds"))

    print(paste0(a/nrow(settings)))
    cat("Time spent = "); cat(difftime(Sys.time(), time.start, units = "mins")); cat("\n")
    cat("Time remain = "); cat(difftime(Sys.time(), time.start, units = "mins")/a*(nrow(settings)-a)); cat("\n")
    return(a.res)

  }
}, mc.cores = ncores)


# summary -----------------------------------------------------------------
num_genes = 10000*10
res = mclapply(1:nrow(settings), function(b){
  a.res = readRDS(paste0(out.dir, "/sim_cate", b, ".rds"))
  a.res$rep = rep(1:10, each = 10000)
  return(a.res)
})

res.tab = do.call(rbind.data.frame, res)
method.cols = 2:8
library(tidyr)
colnames(res.tab) = sapply(colnames(res.tab), function(a){
  switch(a, "aw.type" = "AW", "vda.type" = "VDA", "mod_sel_bic" = "ModSel.BIC", "mod_sel_aic" = "ModSel.AIC",
         "dodr" = "RAIN.DODR", "truth" = "truth", "truth.both" = "truth.both",
         "Sidak_BS" = "Sidak.BS", "Sidak_FS"="Sidak.FS",
         "id" = "id", "rep" = "rep")
})

#split to 10 to calculate sd
res.tab = split(res.tab, res.tab$rep)

class.acc.sep = mclapply(1:10, function(b){
  class.acc = lapply(method.cols, function(a.col){
    a.method.class = as.character(res.tab[[b]][, a.col])
    a.method.class = ifelse(a.method.class%in%c("diffPar", "diffParR2", "diffR2", "same", "change"), "both", a.method.class)
    a.method.tab = table(res.tab[[b]]$truth, a.method.class)
    a.truth.both = gsub("([0-9a-zA-Z])\\..*", "\\1", rownames(a.method.tab))
    a.truth.both = ifelse(a.truth.both%in%c("diffPar", "diffR2", "same"), "both", a.truth.both)
    a.truth.both = ifelse(a.truth.both%in%c("rhy1"), "rhyI", a.truth.both)
    a.method.acc = sapply(1:nrow(a.method.tab), function(a.set){
      a.method.tab[a.set, a.truth.both[a.set]]/sum(a.method.tab[a.set, ])
    })
    out = data.frame(truth = rownames(a.method.tab),
                     truth.both = a.truth.both,
                     acc = a.method.acc)
    colnames(out)[3] = colnames(res.tab[[b]])[a.col]
    return(out)
  })
  class.acc.tab = class.acc[[1]]
  for(i in 2:length(class.acc)){
    class.acc.tab = left_join(class.acc.tab, class.acc[[i]])
  }
  class.acc.tab$rep = b
  return(class.acc.tab)
})
class.acc.all = do.call(rbind.data.frame, class.acc.sep)
write.csv(class.acc.all, paste0(out.dir, "/ClassAccAll.csv"))

library(abind)
class.acc.all2 = class.acc.sep[[1]]
for(i in 2:10){
  class.acc.all2 = abind(class.acc.all2, class.acc.sep[[i]], along = 3)
}
class.acc.mean = cbind(class.acc.sep[[1]][, 1:2],
                       apply(class.acc.all2[, 3:(ncol(class.acc.all2)-1), ], c(1, 2), function(a){mean(as.numeric(a))}))
class.acc.sd = cbind(class.acc.sep[[1]][, 1:2],
                       apply(class.acc.all2[, 3:(ncol(class.acc.all2)-1), ], c(1, 2), function(a){sd(as.numeric(a))}))

# rhyI into both ----------------------------------------------------------
class.both.sep = mclapply(1:10, function(b){
  class.both = lapply(method.cols, function(a.col){
    a.method.class = as.character(res.tab[[b]][, a.col])
    a.method.class = ifelse(a.method.class%in%c("diffPar", "diffParR2", "diffR2", "same", "change"), "both", a.method.class)
    a.method.tab = table(res.tab[[b]]$truth, a.method.class)
    a.truth.both = gsub("([0-9a-zA-Z])\\..*", "\\1", rownames(a.method.tab))
    a.truth.both = ifelse(a.truth.both%in%c("diffPar", "diffR2", "same"), "both", a.truth.both)
    a.truth.both = ifelse(a.truth.both%in%c("rhy1"), "rhyI", a.truth.both)
    a.method.both = sapply(1:nrow(a.method.tab), function(a.set){
      a.method.tab[a.set, "both"]/sum(a.method.tab[a.set, ])
    })
    out = data.frame(truth = rownames(a.method.tab),
                     truth.both = a.truth.both,
                     p.both = a.method.both)
    colnames(out)[3] = colnames(res.tab[[b]])[a.col]
    return(out)
  })
  class.both.tab = class.both[[1]]
  for(i in 2:length(class.both)){
    class.both.tab = left_join(class.both.tab, class.both[[i]])
  }
  class.both.tab$rep = b
  return(class.both.tab)
})
class.both.all = do.call(rbind.data.frame, class.both.sep)
write.csv(class.both.all, paste0(out.dir, "/ClassbothAll.csv"))
class.both.all2 = class.both.sep[[1]]
for(i in 2:10){
  class.both.all2 = abind(class.both.all2, class.both.sep[[i]], along = 3)
}
class.both.mean = cbind(class.both.sep[[1]][, 1:2],
                       apply(class.both.all2[, 3:(ncol(class.both.all2)-1), ], c(1, 2), function(a){mean(as.numeric(a))}))
class.both.sd = cbind(class.both.sep[[1]][, 1:2],
                     apply(class.both.all2[, 3:(ncol(class.both.all2)-1), ], c(1, 2), function(a){sd(as.numeric(a))}))

# rhyI into both P(k=2|k0=1&k>=1)----------------------------------------------------------
class.bothC1.sep = mclapply(1:10, function(b){
  class.bothC1 = lapply(method.cols, function(a.col){
    a.method.class = as.character(res.tab[[b]][, a.col])
    a.method.class = ifelse(a.method.class%in%c("diffPar", "diffParR2", "diffR2", "same", "change"), "both", a.method.class)
    a.method.tab = table(res.tab[[b]]$truth, a.method.class)
    a.truth.both = gsub("([0-9a-zA-Z])\\..*", "\\1", rownames(a.method.tab))
    a.truth.both = ifelse(a.truth.both%in%c("diffPar", "diffR2", "same"), "both", a.truth.both)
    a.truth.both = ifelse(a.truth.both%in%c("rhy1"), "rhyI", a.truth.both)
    a.method.bothC1 = sapply(1:nrow(a.method.tab), function(a.set){
      n.kGTone = sum(a.method.tab[a.set, -which(colnames(a.method.tab)%in%c("arrhy"))])
      a.method.tab[a.set, "both"]/n.kGTone
    })
    out = data.frame(truth = rownames(a.method.tab),
                     truth.both = a.truth.both,
                     bothC1 = a.method.bothC1)
    colnames(out)[3] = colnames(res.tab[[b]])[a.col]
    return(out)
  })
  class.bothC1.tab = class.bothC1[[1]]
  for(i in 2:length(class.bothC1)){
    class.bothC1.tab = left_join(class.bothC1.tab, class.bothC1[[i]])
  }
  class.bothC1.tab$rep = b
  return(class.bothC1.tab)
})
class.bothC1.all = do.call(rbind.data.frame, class.bothC1.sep)
write.csv(class.bothC1.all, paste0(out.dir, "/ClassbothC1All.csv"))
class.bothC1.all2 = class.bothC1.sep[[1]]
for(i in 2:10){
  class.bothC1.all2 = abind(class.bothC1.all2, class.bothC1.sep[[i]], along = 3)
}
class.bothC1.mean = cbind(class.bothC1.sep[[1]][, 1:2],
                       apply(class.bothC1.all2[, 3:(ncol(class.bothC1.all2)-1), ], c(1, 2), function(a){mean(as.numeric(a))}))
class.bothC1.sd = cbind(class.bothC1.sep[[1]][, 1:2],
                     apply(class.bothC1.all2[, 3:(ncol(class.bothC1.all2)-1), ], c(1, 2), function(a){sd(as.numeric(a))}))


##conditional. P(k=2|k0=1&k>=1, M1 = rhyI) #also filter out the cases that M2 = rhyII
class.bothC1M1.sep = mclapply(1:10, function(b){
  class.bothC1M1 = lapply(method.cols, function(a.col){
    a.method.class = as.character(res.tab[[b]][, a.col])
    a.method.class = ifelse(a.method.class%in%c("diffPar", "diffParR2", "diffR2", "same", "change"), "both", a.method.class)
    a.method.tab = table(res.tab[[b]]$truth, a.method.class)
    a.truth.both = gsub("([0-9a-zA-Z])\\..*", "\\1", rownames(a.method.tab))
    a.truth.both = ifelse(a.truth.both%in%c("diffPar", "diffR2", "same"), "both", a.truth.both)
    a.truth.both = ifelse(a.truth.both%in%c("rhy1"), "rhyI", a.truth.both)
    a.method.bothC1M1 = sapply(1:nrow(a.method.tab), function(a.set){
        n.kGTone = sum(a.method.tab[a.set, -which(colnames(a.method.tab)%in%c("rhyII", "arrhy"))])
        a.method.tab[a.set, "both"]/n.kGTone
      })
    out = data.frame(truth = rownames(a.method.tab),
                     truth.both = a.truth.both,
                     bothC1M1 = a.method.bothC1M1)
    colnames(out)[3] = colnames(res.tab[[b]])[a.col]
    return(out)
  })
  class.bothC1M1.tab = class.bothC1M1[[1]]
  for(i in 2:length(class.bothC1M1)){
    class.bothC1M1.tab = left_join(class.bothC1M1.tab, class.bothC1M1[[i]])
  }
  class.bothC1M1.tab$rep = b
  return(class.bothC1M1.tab)
})
class.bothC1M1.all = do.call(rbind.data.frame, class.bothC1M1.sep)
write.csv(class.bothC1M1.all, paste0(out.dir, "/ClassbothC1M1All.csv"))
class.bothC1M1.all2 = class.bothC1M1.sep[[1]]
for(i in 2:10){
  class.bothC1M1.all2 = abind(class.bothC1M1.all2, class.bothC1M1.sep[[i]], along = 3)
}
class.bothC1M1.mean = cbind(class.bothC1M1.sep[[1]][, 1:2],
                          apply(class.bothC1M1.all2[, 3:(ncol(class.bothC1M1.all2)-1), ], c(1, 2), function(a){mean(as.numeric(a))}))
class.bothC1M1.sd = cbind(class.bothC1M1.sep[[1]][, 1:2],
                        apply(class.bothC1M1.all2[, 3:(ncol(class.bothC1M1.all2)-1), ], c(1, 2), function(a){sd(as.numeric(a))}))

# Make plots ------------------------------------------------------
hdir = out.dir
# how to bench mark? #use accuracy to describe and find where is the inaccuracy
# the goal is to classify genes in to rhythmic in I, or both
library(ggplot2)
library(dplyr)
library(tidyr)

## Type I: arrhy to any ----------------------------------------------------
p0 = ggplot(class.acc.mean %>%
              filter(truth=="arrhy") %>%
              pivot_longer(cols = AW:RAIN.DODR, names_to = "Methods", values_to = "acc")%>%
              right_join(class.acc.sd %>%
                           filter(truth=="arrhy") %>%
                           pivot_longer(cols = AW:RAIN.DODR, names_to = "Methods", values_to = "sd"))%>%
              mutate(Methods = factor(Methods, levels = c("AW", "Fisher.BS", "Fisher.FS",
                                                          "Sidak.BS", "Sidak.FS",
                                                          "VDA", "RAIN.DODR", "ModSel.BIC", "ModSel.AIC")))%>%
              filter(!Methods %in% c("Fisher.FS", "Fisher.BS"))
            , aes(x = Methods, y = 1-acc, fill = Methods))+
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=1-acc-sd, ymax=1-acc+sd), width=.2,
                position=position_dodge(.9))+
  geom_hline(yintercept = 0.05, color = "#fcba03", linetype = 2, size = 1, alpha = 0.8)+
  scale_y_continuous(breaks=seq(0, 0.35, by = 0.02))+
  theme_bw()+
  labs(title = expression(A.~TypeI[0]), x = "Methods", y = "Type I error")+
  theme(axis.text.x = element_text(size=10, angle=90))
#= P (Arrhy to rhythmic in one or both)
pdf(paste0(hdir, "/TypeI_error_of_arrhy_to_any_noFisher.pdf"))
print(p0)
dev.off()

## Type I: rhyI to both -----------------------------------------------------
p1 = ggplot(class.both.mean %>%
              filter(grepl("rhy1", truth)) %>%
              pivot_longer(cols = AW:RAIN.DODR, names_to = "Methods", values_to = "acc")%>%
              right_join(class.both.sd %>%
                           filter(grepl("rhy1", truth)) %>%
                           pivot_longer(cols = AW:RAIN.DODR, names_to = "Methods", values_to = "sd"))%>%
              mutate(Methods = factor(Methods, levels = c("AW", "Fisher.BS", "Fisher.FS",
                                                          "Sidak.BS", "Sidak.FS",
                                                          "VDA", "RAIN.DODR", "ModSel.BIC", "ModSel.AIC")))%>%
              mutate(SNR = as.numeric(gsub("rhy1.([\\.0-9]+)", "\\1", truth))) %>%
              filter(!Methods %in% c("VDA", "RAIN.DODR", "ModSel.AIC", "Fisher.FS", "Fisher.BS"))
            , aes(x = SNR, y = acc, color = Methods))+
  geom_line()+
  geom_point(alpha = 0.8, size = 3)+
  geom_errorbar(aes(ymin=acc-sd, ymax=acc+sd), width=.02)+
  geom_hline(yintercept = 0.05, color = "#fcba03", linetype = 2, size = 1, alpha = 0.8)+
  # scale_x_continuous(breaks=seq(0, 1, by = 0.1))+
  scale_y_continuous(breaks=seq(0, 0.25, by = 0.01))+
  theme_bw()+
  labs(title = expression(TypeI[1]), x = expression(A/sigma), y = "Type I error")

pdf(paste0(hdir, "/TypeI_error_of_RhyI_to_Both_noFisher.pdf"))
print(p1)
dev.off()

## Type I cond2: rhyI to both|rhyI -----------------------------------------------------
#The conditional TypeI error: p(k=2|k>=1&k0=1, M1=RhyI)
p1c = ggplot(class.bothC1M1.mean %>%
               filter(grepl("rhy1", truth)) %>%
               pivot_longer(cols = AW:RAIN.DODR, names_to = "Methods", values_to = "acc")%>%
               right_join(class.bothC1M1.sd %>%
                            filter(grepl("rhy1", truth)) %>%
                            pivot_longer(cols = AW:RAIN.DODR, names_to = "Methods", values_to = "sd"))%>%
               mutate(Methods = factor(Methods, levels = c("AW", "Fisher.BS", "Fisher.FS",
                                                           "Sidak.BS", "Sidak.FS",
                                                           "VDA", "RAIN.DODR", "ModSel.BIC", "ModSel.AIC")))%>%
               mutate(SNR = as.numeric(gsub("rhy1.([\\.0-9]+)", "\\1", truth))) %>%
               filter(!Methods %in% c("VDA", "RAIN.DODR", "ModSel.AIC", "Fisher.FS", "Fisher.BS"))
             , aes(x = SNR, y = acc, color = Methods))+
  geom_line()+
  geom_point(alpha = 0.8, size = 3)+
  geom_errorbar(aes(ymin=acc-sd, ymax=acc+sd), width=.02)+
  geom_hline(yintercept = 0.05, color = "#fcba03", linetype = 2, size = 1, alpha = 0.8)+
  # scale_x_continuous(breaks=seq(0, 1, by = 0.1))+
  scale_y_continuous(breaks=seq(0, 0.25, by = 0.01))+
  theme_bw()+
  labs(title = expression(B.~cTypeI[1]), x = expression(A/sigma), y = "Type I error")

pdf(paste0(hdir, "/TypeI_error_of_RhyI_to_BothC1M1_noFisher.pdf.pdf"))
print(p1c)
dev.off()

pdf(paste0(hdir, "/TypeI_error_noFisher.pdf"), width = 10, height = 5)
gridExtra::grid.arrange(p0, p1c, ncol = 2)
dev.off()

## power2 -------------------------------------------------------------------
p1.1 = ggplot(class.acc.mean %>%
                filter(grepl("rhy1", truth)) %>%
                pivot_longer(cols = AW:RAIN.DODR, names_to = "Methods", values_to = "acc")%>%
                right_join(class.acc.sd %>%
                             filter(grepl("rhy1", truth)) %>%
                             pivot_longer(cols = AW:RAIN.DODR, names_to = "Methods", values_to = "sd"))%>%
                mutate(Methods = factor(Methods, levels = c("AW", "Fisher.BS", "Fisher.FS",
                                                            "Sidak.BS", "Sidak.FS",
                                                            "VDA", "RAIN.DODR", "ModSel.BIC", "ModSel.AIC")))%>%
                mutate(SNR = as.numeric(gsub("rhy1.([\\.0-9]+)", "\\1", truth))) %>%
                filter(!Methods %in% c("AW","ModSel.AIC", "RAIN.DODR", "VDA", "Fisher.FS", "ModSel.BIC", "Fisher.BS"))
              , aes(x = SNR, y = acc, color = Methods))+
  geom_line()+
  geom_point(alpha = 0.8, size = 3)+
  scale_x_continuous(breaks=seq(0, 1, by = 0.1))+
  scale_y_continuous(breaks=seq(0.2, 1, by = 0.1))+
  # coord_cartesian(xlim=c(0.4, 1), ylim=c(0.15, 0.90))+
  theme_bw()+
  labs(title = expression(A.~Power[1]), x = expression(A/sigma), y = "Power")
# the inaccuracy went to arrhy because of the low power.

p2.same = ggplot(class.acc.mean %>%
                   filter(grepl("same", truth)) %>%
                   pivot_longer(cols = AW:RAIN.DODR, names_to = "Methods", values_to = "acc")%>%
                   right_join(class.acc.sd %>%
                                filter(grepl("same", truth)) %>%
                                pivot_longer(cols = AW:RAIN.DODR, names_to = "Methods", values_to = "sd"))%>%
                   mutate(Methods = factor(Methods, levels = c("AW", "Fisher.BS", "Fisher.FS",
                                                               "Sidak.BS", "Sidak.FS",
                                                               "VDA", "RAIN.DODR", "ModSel.BIC", "ModSel.AIC")))%>%
                   mutate(SNR = as.numeric(gsub("same.([\\.0-9]+)", "\\1", truth))) %>%
                   filter(!Methods %in% c("AW","ModSel.AIC", "RAIN.DODR", "VDA", "Fisher.FS", "ModSel.BIC", "Fisher.BS"))
                 , aes(x = SNR, y = acc, color = Methods))+
  geom_line()+
  geom_point(alpha = 0.8, size = 3)+
  scale_x_continuous(breaks=seq(0, 1, by = 0.1))+
  scale_y_continuous(breaks=seq(0.2, 1, by = 0.1))+
  # coord_cartesian(xlim=c(0.4, 1))+
  theme_bw()+
  labs(title = expression(B.~Power[2]~(Same~SNR)), x = expression(A/sigma), y = "Power")

p2.diffPar = ggplot(class.acc.mean %>%
                      filter(grepl("diffPar", truth)) %>%
                      pivot_longer(cols = AW:RAIN.DODR, names_to = "Methods", values_to = "acc")%>%
                      right_join(class.acc.sd %>%
                                   filter(grepl("diffPar", truth)) %>%
                                   pivot_longer(cols = AW:RAIN.DODR, names_to = "Methods", values_to = "sd"))%>%
                      mutate(Methods = factor(Methods, levels = c("AW", "Fisher.BS", "Fisher.FS",
                                                                  "Sidak.BS", "Sidak.FS",
                                                                  "VDA", "RAIN.DODR", "ModSel.BIC", "ModSel.AIC")))%>%
                      mutate(SNR = as.numeric(gsub("diffPar.([\\.0-9]+)", "\\1", truth))) %>%
                      filter(!Methods %in% c("AW","ModSel.AIC", "RAIN.DODR", "VDA", "Fisher.FS", "ModSel.BIC", "Fisher.BS"))
                    , aes(x = SNR, y = acc, color = Methods))+
  geom_line()+
  geom_point(alpha = 0.8, size = 3)+
  # scale_x_continuous(breaks=c(pi/8, pi/4, pi/2, pi), labels = c("pi/8", "pi/4", "pi/2", "pi"))+
  scale_y_continuous(breaks=seq(0.75, 1, by = 0.01))+
  theme_bw()+
  labs(title = expression(C.~Power[2]~(Same~SNR~with~different~phase)), x = expression(Delta~phi), y = "Power")

p2.diffR2 = ggplot(class.acc.mean %>%
                     filter(grepl("diffR", truth)) %>%
                     pivot_longer(cols = AW:RAIN.DODR, names_to = "Methods", values_to = "acc")%>%
                     right_join(class.acc.sd %>%
                                  filter(grepl("diffR", truth)) %>%
                                  pivot_longer(cols = AW:RAIN.DODR, names_to = "Methods", values_to = "sd"))%>%
                     mutate(Methods = factor(Methods, levels = c("AW", "Fisher.BS", "Fisher.FS",
                                                                 "Sidak.BS", "Sidak.FS",
                                                                 "VDA", "RAIN.DODR", "ModSel.BIC", "ModSel.AIC")))%>%
                     mutate(SNR = 1/as.numeric(gsub("diffR2.([\\.0-9]+)", "\\1", truth))) %>%
                     filter(!Methods %in% c("AW","ModSel.AIC", "RAIN.DODR", "VDA", "Fisher.FS", "ModSel.BIC", "Fisher.BS"))
                   , aes(x = SNR, y = acc, color = Methods))+
  geom_line()+
  geom_point(alpha = 0.8, size = 3)+
  # scale_x_continuous(breaks=c(0.67, 0.83, 1.25, 1.67, 2.50), labels = c("0.67", "0.83", "1.25", "1.67", "2.50"))+
  scale_x_continuous(breaks=c(1.50, 1.2, 0.8, 0.6, 0.4), labels = c("1.50", "1.2", "0.8", "0.6", "0.4"))+
  scale_y_continuous(breaks=seq(0, 1, by = 0.1))+
  theme_bw()+
  labs(title = expression(D.~Power[2]~(Different~SNR)), x = expression(A[2]/sigma[2]), y = "Power")

pdf(paste0(hdir, "/Power_noFisher.pdf"), width = 10)
gridExtra::grid.arrange(p1.1, p2.same, p2.diffPar, p2.diffR2, ncol = 2)
dev.off()

