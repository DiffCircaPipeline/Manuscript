rm(list = ls())
hdir = "/home/xix66/circadian/ThePipeline/CircadianPipeline_private/LR_overall" #speed3
library(parallel)
out.dir = paste0(hdir, "/output_DRparam")
dir.create(out.dir, recursive = TRUE)


# A&phase&M: type I -------------------------------------------
#simulation setting
num_genes = 10000*10
n.sample = 30
A1 = 3
phase1 = 0
M1 = 5
sigma1 = 1

settings = expand.grid(A = c(2, 3, 4),
                       phase = c(0, 6, 12),
                       M = c(5, 6, 7),
                       sigma = 1)
settings.keep = sapply(1:nrow(settings), function(i){
  a.settings = settings[i, ]
  return( sum(a.settings$A==3, a.settings$phase ==0, a.settings$M==5)>=3)
})
settings = settings[settings.keep, ]

time.start = Sys.time()
res = lapply(1:nrow(settings), function(a){
  a.sigma = settings$sigma[a]
  a.M = settings$M[a]
  a.phase = settings$phase[a]
  a.A = settings$A[a]

  B.res = mclapply(1:num_genes, function(b){
    set.seed(b)
    # sample times
    x1.time = runif(n.sample, min=0, max=24)
    x2.time = runif(n.sample, min=0, max=24)

    # noise
    noise.mat1 = matrix(rnorm(1*n.sample, 0, a.sigma), ncol=n.sample, nrow=1)
    noise.mat2 = matrix(rnorm(1*n.sample, 0, a.sigma), ncol=n.sample, nrow=1)

    # signal
    signal.mat1 = t(sapply(1:1, function(g){a.M + a.A * cos((2*pi/24)*(x1.time - a.phase))}))
    signal.mat2 = t(sapply(1:1, function(g){a.M + a.A * cos((2*pi/24)*(x2.time - a.phase))}))

    t.all = c(x1.time, x2.time)
    y.all = c(noise.mat1 + signal.mat1, noise.mat2 + signal.mat2)
    group.all = c(rep(0, length(x1.time)), rep(1, length(x2.time)))
    x1 = list(data = data.frame(noise.mat1 + signal.mat1),
              time = x1.time,
              gname = paste("gene", seq_len(nrow(signal.mat1))))
    x2 = list(data = data.frame(noise.mat2 + signal.mat2),
              time = x2.time,
              gname = paste("gene", seq_len(nrow(signal.mat2))))
    DCP.obj = DiffCircaPipeline::DCP_Rhythmicity(x1, x2)
    DCP.Dparam = DiffCircaPipeline::DCP_DiffPar(DCP.obj, "A&phase&M", "both")
    testA = diffCircadian::LR_diff(x1.time, noise.mat1 + signal.mat1, x2.time, noise.mat2 + signal.mat2, period = 24, FN = TRUE, type="amplitude")
    testphase = diffCircadian::LR_diff(x1.time, noise.mat1 + signal.mat1, x2.time, noise.mat2 + signal.mat2, period = 24, FN = TRUE, type="phase")
    testM = diffCircadian::LR_diff(x1.time, noise.mat1 + signal.mat1, x2.time, noise.mat2 + signal.mat2, period = 24, FN = TRUE, type="basal")

    a.tab.sidak = data.frame(Combined = DCP.Dparam$p.overall<0.05,
                             A.w.overall = testA$pvalue<0.05, phase.w.overall = testphase$pvalue<0.05,
                             M.w.overall = testM$pvalue<0.05,
                             any.w.P = testA$pvalue<0.05|testphase$pvalue<0.05|testM$pvalue<0.05,
                             A.o.overall = DCP.Dparam$post.hoc.A.By.p,
                             phase.o.overall = DCP.Dparam$post.hoc.peak.By.p,
                             M.o.overall = DCP.Dparam$post.hoc.M.By.p,
                             any.o.overall = DCP.Dparam$post.hoc.A.By.p|DCP.Dparam$post.hoc.peak.By.p|DCP.Dparam$post.hoc.M.By.p,
                             setting = a,
                             method = "Sidak",
                             sigma = settings$sigma[a],
                             M = settings$M[a],
                             phase = settings$phase[a],
                             A = settings$A[a],
                             test.ind = b)

    return(a.tab.sidak)

  }, mc.cores = 40)

  diff.res = do.call(rbind.data.frame, B.res)

  save(diff.res, file=paste0(out.dir,  "/LR_overall_parameterAPM_simTypeI", "_A", a.A, "_m", a.M, "_phase", round(a.phase, 2), "_sigma", a.sigma,  ".RData"))
  #mark time
  print(paste0(a/nrow(settings)))
  cat("Time spent = "); cat(difftime(Sys.time(), time.start, units = "mins")); cat("\n")
  cat("Time remain = "); cat(difftime(Sys.time(), time.start, units = "mins")/a*(nrow(settings)-a)); cat("\n")

})

# A&P: type I ---------------------------------------------------------
rm(list = ls())
hdir = "/home/xix66/circadian/ThePipeline/CircadianPipeline_private/LR_overall" #speed3
out.dir = paste0(hdir, "/output_DRparam")

num_genes = 10000*10
n.sample = 30
A1 = 3
phase1 = 0
M1 = 5
sigma1 = 1

settings = expand.grid(A = c(2, 3, 4),
                       phase = c(0, 6, 12),
                       M = c(5, 6, 7),
                       sigma = 1)
settings.keep = sapply(1:nrow(settings), function(i){
  a.settings = settings[i, ]
  return( sum(a.settings$A==3, a.settings$phase ==0, a.settings$M==5)>=3)
})
settings = settings[settings.keep, ]

time.start = Sys.time()
res = lapply(1:nrow(settings), function(a){
  a.sigma = settings$sigma[a]
  a.M = settings$M[a]
  a.phase = settings$phase[a]
  a.A = settings$A[a]

  B.res = mclapply(1:num_genes, function(b){
    set.seed(b)
    # sample times
    x1.time = runif(n.sample, min=0, max=24)
    x2.time = runif(n.sample, min=0, max=24)

    # noise
    noise.mat1 = matrix(rnorm(1*n.sample, 0, a.sigma), ncol=n.sample, nrow=1)
    noise.mat2 = matrix(rnorm(1*n.sample, 0, a.sigma), ncol=n.sample, nrow=1)

    # signal
    signal.mat1 = t(sapply(1:1, function(g){a.M + a.A * cos((2*pi/24)*(x1.time - a.phase))}))
    signal.mat2 = t(sapply(1:1, function(g){a.M + a.A * cos((2*pi/24)*(x2.time - a.phase))}))

    t.all = c(x1.time, x2.time)
    y.all = c(noise.mat1 + signal.mat1, noise.mat2 + signal.mat2)
    group.all = c(rep(0, length(x1.time)), rep(1, length(x2.time)))
    x1 = list(data = data.frame(noise.mat1 + signal.mat1),
              time = x1.time,
              gname = paste("gene", seq_len(nrow(signal.mat1))))
    x2 = list(data = data.frame(noise.mat2 + signal.mat2),
              time = x2.time,
              gname = paste("gene", seq_len(nrow(signal.mat2))))
    DCP.obj = DiffCircaPipeline::DCP_Rhythmicity(x1, x2)
    DCP.Dparam = DiffCircaPipeline::DCP_DiffPar(DCP.obj, "A&phase", "both")
    testA = diffCircadian::LR_diff(x1.time, noise.mat1 + signal.mat1, x2.time, noise.mat2 + signal.mat2, period = 24, FN = TRUE, type="amplitude")
    testphase = diffCircadian::LR_diff(x1.time, noise.mat1 + signal.mat1, x2.time, noise.mat2 + signal.mat2, period = 24, FN = TRUE, type="phase")

    a.tab.sidak = data.frame(Combined = DCP.Dparam$p.overall<0.05,
                             A.w.overall = testA$pvalue<0.05, phase.w.overall = testphase$pvalue<0.05,
                             any.w.P = testA$pvalue<0.05|testphase$pvalue<0.05,
                             A.o.overall = DCP.Dparam$post.hoc.A.By.p,
                             phase.o.overall = DCP.Dparam$post.hoc.peak.By.p,
                             any.o.overall = DCP.Dparam$post.hoc.A.By.p|DCP.Dparam$post.hoc.peak.By.p,
                             setting = a,
                             method = "Sidak",
                             sigma = settings$sigma[a],
                             M = settings$M[a],
                             phase = settings$phase[a],
                             A = settings$A[a],
                             test.ind = b)

    return(a.tab.sidak)

  }, mc.cores = 40)

  diff.res = do.call(rbind.data.frame, B.res)

  save(diff.res, file=paste0(out.dir,  "/LR_overall_parameterAP_simTypeI", "_A", a.A, "_m", a.M, "_phase", round(a.phase, 2), "_sigma", a.sigma,  ".RData"))
  #mark time
  print(paste0(a/nrow(settings)))
  cat("Time spent = "); cat(difftime(Sys.time(), time.start, units = "mins")); cat("\n")
  cat("Time remain = "); cat(difftime(Sys.time(), time.start, units = "mins")/a*(nrow(settings)-a)); cat("\n")

})

# Summary -----------------------------------------------------------------
library(dplyr)
library(tidyr)
settings = expand.grid(A = c(2, 3, 4),
                       phase = c(0, 6, 12),
                       M = c(5, 6, 7),
                       sigma = 1)
settings.keep = sapply(1:nrow(settings), function(i){
  a.settings = settings[i, ]
  return( sum(a.settings$A==3, a.settings$phase ==0, a.settings$M==5)>=3)
})
settings = settings[settings.keep, ]

res = mclapply(1:nrow(settings), function(a){
  a.sigma = settings$sigma[a]
  a.M = settings$M[a]
  a.phase = settings$phase[a]
  a.A = settings$A[a]

  load(file=paste0(out.dir,  "/LR_overall_parameterAPM_simTypeI", "_A", a.A, "_m", a.M, "_phase", round(a.phase, 2), "_sigma", a.sigma,  ".RData"))
  return(diff.res)
  #mark time
},mc.cores = 10)

rep = rep(1:10, each = 10000)
res = lapply(res, function(a.res){
  a.res$rep = rep
  return(a.res)
})

res.tab = do.call(rbind.data.frame, res)
res.tab.mean.s = res.tab %>%
  split(f = list(res.tab$M, res.tab$phase, res.tab$A, res.tab$rep, res.tab$method)) %>%
  lapply(function(x){
    a.tab = data.frame(Combined = mean(x$Combined),
                       any.sep = mean(x$any.w.P),
                       A.o.overall = mean(x$A.o.overall), phase.o.overall = mean(x$phase.o.overall), M.o.overall = mean(x$M.o.overall),
                       any.o.overall = mean(x$any.o.overall),
                       method = x$method[1],
                       sigma = x$sigma[1],
                       M = x$M[1],
                       phase = x$phase[1],
                       A = x$A[1])
    return(a.tab)
  })
res.tab.mean = do.call(rbind.data.frame, res.tab.mean.s)
res.tab.mean = res.tab.mean[complete.cases(res.tab.mean), ]
res.tab.mean$M = factor(res.tab.mean$M)
res.tab.mean$A = factor(res.tab.mean$A)
res.tab.mean$phase = factor(round(res.tab.mean$phase,2))
res.tab.mean$sigma = factor(round(res.tab.mean$sigma,2))

res.tab.mean = res.tab.mean %>% pivot_longer(cols = c(1:6),
                                             names_to = "Tests",
                                             values_to = "Type I error rate")


## plot version 1 ----------------------------------------------------------
res.tab.mean$Tests = sapply(res.tab.mean$Tests, function(a){
  switch(a,
         "Combined" = "Combined",
         "A.o.overall" = "Post hoc A",
         "phase.o.overall" = "Post hoc phase",
         "M.o.overall" = "Post hoc M",
         "any.o.overall" = "Post hoc any",
         "any.sep" = "Seperately any")
})
res.tab.mean$Tests = factor(res.tab.mean$Tests, levels = c("Combined", "Post hoc A", "Post hoc phase", "Post hoc M", "Post hoc any", "Seperately any"))

library(ggplot2)
res.tab.mean2 = res.tab.mean[, -c(1)] %>%
  split(f = list(res.tab.mean$Tests)) %>%
  lapply(function(x){
    x.tab = data.frame(sigma = x$sigma[1],
                       M = x$M[1],
                       phase = x$phase[1],
                       A = x$A[1],
                       mean = mean(x$`Type I error rate`),
                       sd = sd(x$`Type I error rate`),
                       Tests = x$Tests[1])
    return(x.tab)
  })
res.tab.mean2 = do.call(rbind.data.frame, res.tab.mean2)

pTypeI3 = ggplot(res.tab.mean2, aes(x = Tests, y = mean, fill = Tests))+
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  geom_hline(yintercept = 0.05, color = "#fcba03", linetype = 2, size = 1, alpha = 0.8)+
  # scale_y_continuous(breaks=seq(0, 0.35, by = 0.02))+
  theme_bw()+
  scale_x_discrete(labels = c(expression(Combined~(Delta*A*","~Delta*phi*","~Delta*M)),
                              expression(Post*" "*hoc~Delta*A),
                              expression(Post*" "*hoc~Delta*phi),
                              expression(Post*" "*hoc~Delta*M),
                              expression(Post*" "*hoc~(Delta*A*","~Delta*phi*","~Delta*M)),
                              expression(Separate~(Delta*A*","~Delta*phi*","~Delta*M))))+
  scale_fill_manual(labels = c(expression(Combined~(Delta*A*","~Delta*phi*","~Delta*M)),
                               expression(Post*" "*hoc~Delta*A),
                               expression(Post*" "*hoc~Delta*phi),
                               expression(Post*" "*hoc~Delta*M),
                               expression(Post*" "*hoc~(Delta*A*","~Delta*phi*","~Delta*M)),
                               expression(Separate~(Delta*A*","~Delta*phi*","~Delta*M))),
                    # values = c("#F8766D", "#D39200", "#00BA38", "#00B9E3", "#DB72FB"))+
                    values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"))+
  labs(title = expression(B*.~Type~I~error~(H[0][","]["("][A][", "][phi][", "][M][")"])),
       subtitle = expression((A*"="*3*","~phi*"="*0*","~M*"="*5)),
       x = "Tests", y = "Type I error")+
  theme(axis.text.x = element_text(angle=90), legend.text.align = 0)

# parameters A&phase
settings = expand.grid(A = c(2, 3, 4),
                       phase = c(0, 6, 12),
                       M = c(5, 6, 7),
                       sigma = 1)
settings.keep = sapply(1:nrow(settings), function(i){
  a.settings = settings[i, ]
  return( sum(a.settings$A==3, a.settings$phase ==0, a.settings$M==5)>=3)
})
settings = settings[settings.keep, ]

res = mclapply(1:nrow(settings), function(a){
  a.sigma = settings$sigma[a]
  a.M = settings$M[a]
  a.phase = settings$phase[a]
  a.A = settings$A[a]

  load(file=paste0(out.dir,  "/LR_overall_parameterAP_simTypeI", "_A", a.A, "_m", a.M, "_phase", round(a.phase, 2), "_sigma", a.sigma,  ".RData"))
  return(diff.res)
  #mark time
},mc.cores = 10)

rep = rep(1:10, each = 10000)
res = lapply(res, function(a.res){
  a.res$rep = rep
  return(a.res)
})

res.tab = do.call(rbind.data.frame, res)
res.tab.mean.s = res.tab %>%
  split(f = list(res.tab$M, res.tab$phase, res.tab$A, res.tab$rep, res.tab$method)) %>%
  lapply(function(x){
    a.tab = data.frame(Combined = mean(x$Combined),
                       any.sep = mean(x$any.w.P),
                       A.o.overall = mean(x$A.o.overall), phase.o.overall = mean(x$phase.o.overall),
                       any.o.overall = mean(x$any.o.overall),
                       method = x$method[1],
                       sigma = x$sigma[1],
                       M = x$M[1],
                       phase = x$phase[1],
                       A = x$A[1])
    return(a.tab)
  })
res.tab.mean = do.call(rbind.data.frame, res.tab.mean.s)
res.tab.mean = res.tab.mean[complete.cases(res.tab.mean), ]
res.tab.mean$M = factor(res.tab.mean$M)
res.tab.mean$A = factor(res.tab.mean$A)
res.tab.mean$phase = factor(round(res.tab.mean$phase,2))
res.tab.mean$sigma = factor(round(res.tab.mean$sigma,2))

res.tab.mean = res.tab.mean %>% pivot_longer(cols = c(1:5),
                                             names_to = "Tests",
                                             values_to = "Type I error rate")
res.tab.mean$Tests = sapply(res.tab.mean$Tests, function(a){
  switch(a,
         "Combined" = "Combined",
         "A.o.overall" = "Post hoc A",
         "phase.o.overall" = "Post hoc phase",
         "any.o.overall" = "Post hoc any",
         "any.sep" = "Seperately any")
})
res.tab.mean$Tests = factor(res.tab.mean$Tests, levels = c("Combined", "Post hoc A", "Post hoc phase", "Post hoc any", "Seperately any"))

library(ggplot2)
res.tab.mean2 = res.tab.mean[, -c(1)] %>%
  split(f = list(res.tab.mean$Tests)) %>%
  lapply(function(x){
    x.tab = data.frame(sigma = x$sigma[1],
                       M = x$M[1],
                       phase = x$phase[1],
                       A = x$A[1],
                       mean = mean(x$`Type I error rate`),
                       sd = sd(x$`Type I error rate`),
                       Tests = x$Tests[1])
    return(x.tab)
  })
res.tab.mean2 = do.call(rbind.data.frame, res.tab.mean2)

pTypeI2 = ggplot(res.tab.mean2, aes(x = Tests, y = mean, fill = Tests))+
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  geom_hline(yintercept = 0.05, color = "#fcba03", linetype = 2, size = 1, alpha = 0.8)+
  # scale_y_continuous(breaks=seq(0, 0.35, by = 0.02))+
  theme_bw()+
  scale_x_discrete(labels = c(expression(Combined~(Delta*A*","~Delta*phi)),
                              expression(Post*" "*hoc~Delta*A),
                              expression(Post*" "*hoc~Delta*phi),
                              expression(Post*" "*hoc~(Delta*A*","~Delta*phi)),
                              expression(Separate~(Delta*A*","~Delta*phi))))+
  scale_fill_manual(labels = c(expression(Combined~(Delta*A*","~Delta*phi)),
                               expression(Post*" "*hoc~Delta*A),
                               expression(Post*" "*hoc~Delta*phi),
                               expression(Post*" "*hoc~(Delta*A*","~Delta*phi)),
                               expression(Separate~(Delta*A*","~Delta*phi))),
                     values = c("#F8766D", "#D39200", "#00BA38", "#00B9E3", "#DB72FB"))+
                    # values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"))+
  labs(title = expression(A*.~Type~I~error~(H[0][","]["("][A][", "][phi][")"])),
       subtitle = expression((A*"="*3*","~phi*"="*0)),
       x = "Tests", y = "Type I error")+
  theme(axis.text.x = element_text(angle=90), legend.text.align = 0)

pdf(paste0(out.dir, "/paper_LR_overall_TypeI.pdf"), width = 10, height = 5)
lay <- matrix(c(rep(1, 19), rep(3, 1), rep(2, 20)), ncol = 2, byrow = FALSE)
gridExtra::grid.arrange(grobs = list(pTypeI2, pTypeI3), layout_matrix = lay)
dev.off()

## plot version 2 ----------------------------------------------------------
library(dplyr)
library(tidyr)
settings = expand.grid(A = c(2, 3, 4),
                       phase = c(0, 6, 12),
                       M = c(5, 6, 7),
                       sigma = 1)
settings.keep = sapply(1:nrow(settings), function(i){
  a.settings = settings[i, ]
  return( sum(a.settings$A==3, a.settings$phase ==0, a.settings$M==5)>=3)
})
settings = settings[settings.keep, ]

res = mclapply(1:nrow(settings), function(a){
  a.sigma = settings$sigma[a]
  a.M = settings$M[a]
  a.phase = settings$phase[a]
  a.A = settings$A[a]

  load(file=paste0(out.dir,  "/LR_overall_parameterAPM_simTypeI", "_A", a.A, "_m", a.M, "_phase", round(a.phase, 2), "_sigma", a.sigma,  ".RData"))
  return(diff.res)
  #mark time
},mc.cores = 10)

rep = rep(1:10, each = 10000)
res = lapply(res, function(a.res){
  a.res$rep = rep
  return(a.res)
})

res.tab = do.call(rbind.data.frame, res)
res.tab.mean.s = res.tab %>%
  split(f = list(res.tab$M, res.tab$phase, res.tab$A, res.tab$rep, res.tab$method)) %>%
  lapply(function(x){
    a.tab = data.frame(Combined = mean(x$Combined),
                       any.sep = mean(x$any.w.P),
                       A.o.overall = mean(x$A.o.overall), phase.o.overall = mean(x$phase.o.overall), M.o.overall = mean(x$M.o.overall),
                       any.o.overall = mean(x$any.o.overall),
                       method = x$method[1],
                       sigma = x$sigma[1],
                       M = x$M[1],
                       phase = x$phase[1],
                       A = x$A[1])
    return(a.tab)
  })
res.tab.mean = do.call(rbind.data.frame, res.tab.mean.s)
res.tab.mean = res.tab.mean[complete.cases(res.tab.mean), ]
res.tab.mean$M = factor(res.tab.mean$M)
res.tab.mean$A = factor(res.tab.mean$A)
res.tab.mean$phase = factor(round(res.tab.mean$phase,2))
res.tab.mean$sigma = factor(round(res.tab.mean$sigma,2))

res.tab.mean = res.tab.mean %>% pivot_longer(cols = c(1:6),
                                             names_to = "Tests",
                                             values_to = "Type I error rate")

res.tab.mean$Tests = sapply(res.tab.mean$Tests, function(a){
  switch(a,
         "Combined" = "Combined",
         "A.o.overall" = "Post hoc A",
         "phase.o.overall" = "Post hoc phase",
         "M.o.overall" = "Post hoc M",
         "any.o.overall" = "Post hoc any",
         "any.sep" = "Seperately any")
})
res.tab.mean$Tests = factor(res.tab.mean$Tests, levels = c("Combined", "Post hoc A", "Post hoc phase", "Post hoc M", "Post hoc any", "Seperately any"))
#          "Combined" = "Fisrt-stage combined test",
#          "A.o.overall" = "Second-stage post hoc A",
#          "phase.o.overall" = "Second-stage post hoc phase",
#          "M.o.overall" = "Second-stage post hoc M",
#          "any.o.overall" = "Union of second-stage tests",
#          "any.sep" = "Seperate tests without multiplicity correction")
library(ggplot2)
res.tab.mean2 = res.tab.mean[, -c(1)] %>%
  split(f = list(res.tab.mean$Tests)) %>%
  lapply(function(x){
    x.tab = data.frame(sigma = x$sigma[1],
                       M = x$M[1],
                       phase = x$phase[1],
                       A = x$A[1],
                       mean = mean(x$`Type I error rate`),
                       sd = sd(x$`Type I error rate`),
                       Tests = x$Tests[1])
    return(x.tab)
  })
res.tab.mean2 = do.call(rbind.data.frame, res.tab.mean2)

pTypeI3 = ggplot(res.tab.mean2, aes(x = Tests, y = mean, fill = Tests))+
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  geom_hline(yintercept = 0.05, color = "#fcba03", linetype = 2, size = 1, alpha = 0.8)+
  # scale_y_continuous(breaks=seq(0, 0.35, by = 0.02))+
  theme_bw()+
  # scale_x_discrete(labels = c(expression(First*"-"*stage~combined~test),
  #                             expression(Second*"-"*stage~post*" "*hoc~Delta*A),
  #                             expression(Second*"-"*stage~post*" "*hoc~Delta*phi),
  #                             expression(Second*"-"*stage~post*" "*hoc~Delta*M),
  #                             expression(Union~of~"second"*"-"*stage~tests),
  #                             expression(Seperate~tests~without~multiplicity~correction)))+
  scale_fill_manual(labels = c(expression(First*"-"*stage~global~test),
                               expression(Second*"-"*stage~post*" "*hoc~Delta*A),
                               expression(Second*"-"*stage~post*" "*hoc~Delta*phi),
                               expression(Second*"-"*stage~post*" "*hoc~Delta*M),
                               expression(Union~of~"second"*"-"*stage~tests),
                               expression(Seperate~tests~without~multiplicity~correction)),
                    # values = c("#F8766D", "#D39200", "#00BA38", "#00B9E3", "#DB72FB"))+
                    values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"))+
  labs(title = expression(B*.~Type~I~error~(H[0][","]["("][A][", "][phi][", "][M][")"])),
       subtitle = expression((A*"="*3*","~phi*"="*0*","~M*"="*5)),
       x = "Tests", y = "Type I error")+
  theme(axis.text.x = element_blank(),axis.ticks = element_blank(), legend.text.align = 0)

# parameters A&phase
settings = expand.grid(A = c(2, 3, 4),
                       phase = c(0, 6, 12),
                       M = c(5, 6, 7),
                       sigma = 1)
settings.keep = sapply(1:nrow(settings), function(i){
  a.settings = settings[i, ]
  return( sum(a.settings$A==3, a.settings$phase ==0, a.settings$M==5)>=3)
})
settings = settings[settings.keep, ]

res = mclapply(1:nrow(settings), function(a){
  a.sigma = settings$sigma[a]
  a.M = settings$M[a]
  a.phase = settings$phase[a]
  a.A = settings$A[a]

  load(file=paste0(out.dir,  "/LR_overall_parameterAP_simTypeI", "_A", a.A, "_m", a.M, "_phase", round(a.phase, 2), "_sigma", a.sigma,  ".RData"))
  return(diff.res)
  #mark time
},mc.cores = 10)

rep = rep(1:10, each = 10000)
res = lapply(res, function(a.res){
  a.res$rep = rep
  return(a.res)
})

res.tab = do.call(rbind.data.frame, res)
res.tab.mean.s = res.tab %>%
  split(f = list(res.tab$M, res.tab$phase, res.tab$A, res.tab$rep, res.tab$method)) %>%
  lapply(function(x){
    a.tab = data.frame(Combined = mean(x$Combined),
                       any.sep = mean(x$any.w.P),
                       A.o.overall = mean(x$A.o.overall), phase.o.overall = mean(x$phase.o.overall),
                       any.o.overall = mean(x$any.o.overall),
                       method = x$method[1],
                       sigma = x$sigma[1],
                       M = x$M[1],
                       phase = x$phase[1],
                       A = x$A[1])
    return(a.tab)
  })
res.tab.mean = do.call(rbind.data.frame, res.tab.mean.s)
res.tab.mean = res.tab.mean[complete.cases(res.tab.mean), ]
res.tab.mean$M = factor(res.tab.mean$M)
res.tab.mean$A = factor(res.tab.mean$A)
res.tab.mean$phase = factor(round(res.tab.mean$phase,2))
res.tab.mean$sigma = factor(round(res.tab.mean$sigma,2))

res.tab.mean = res.tab.mean %>% pivot_longer(cols = c(1:5),
                                             names_to = "Tests",
                                             values_to = "Type I error rate")
res.tab.mean$Tests = sapply(res.tab.mean$Tests, function(a){
  switch(a,
         "Combined" = "Combined",
         "A.o.overall" = "Post hoc A",
         "phase.o.overall" = "Post hoc phase",
         "any.o.overall" = "Post hoc any",
         "any.sep" = "Seperately any")
})
res.tab.mean$Tests = factor(res.tab.mean$Tests, levels = c("Combined", "Post hoc A", "Post hoc phase", "Post hoc any", "Seperately any"))

library(ggplot2)
res.tab.mean2 = res.tab.mean[, -c(1)] %>%
  split(f = list(res.tab.mean$Tests)) %>%
  lapply(function(x){
    x.tab = data.frame(sigma = x$sigma[1],
                       M = x$M[1],
                       phase = x$phase[1],
                       A = x$A[1],
                       mean = mean(x$`Type I error rate`),
                       sd = sd(x$`Type I error rate`),
                       Tests = x$Tests[1])
    return(x.tab)
  })
res.tab.mean2 = do.call(rbind.data.frame, res.tab.mean2)

pTypeI2 = ggplot(res.tab.mean2, aes(x = Tests, y = mean, fill = Tests))+
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  geom_hline(yintercept = 0.05, color = "#fcba03", linetype = 2, size = 1, alpha = 0.8)+
  # scale_y_continuous(breaks=seq(0, 0.35, by = 0.02))+
  theme_bw()+
  # scale_x_discrete(labels = c(expression(First*"-"*stage~combined~test),
  #                             expression(Second*"-"*stage~post*" "*hoc~Delta*A),
  #                             expression(Second*"-"*stage~post*" "*hoc~Delta*phi),
  #                             expression(Union~of~"second"*"-"*stage*tests),
  #                             expression(Seperate~tests~without~multiplicity~correction)))+
  scale_fill_manual(labels = c(expression(First*"-"*stage~global~test),
                               expression(Second*"-"*stage~post*" "*hoc~Delta*A),
                               expression(Second*"-"*stage~post*" "*hoc~Delta*phi),
                               expression(Union~of~"second"*"-"*stage~tests),
                               expression(Seperate~tests~without~multiplicity~correction)),
                    values = c("#F8766D", "#D39200", "#00BA38", "#619CFF", "#F564E3"))+
  # values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"))+
  labs(title = expression(A*.~Type~I~error~(H[0][","]["("][A][", "][phi][")"])),
       subtitle = expression((A*"="*3*","~phi*"="*0)),
       x = "Tests", y = "Type I error")+
  theme(axis.text.x = element_blank(),axis.ticks = element_blank(), legend.text.align = 0)

pdf(paste0(out.dir, "/paper_LR_overall_TypeI_renamed.pdf"), width = 12, height = 3.5)
lay <- matrix(c(rep(1, 19), rep(3, 1), rep(2, 20)), ncol = 2, byrow = FALSE)
gridExtra::grid.arrange(grobs = list(pTypeI2, pTypeI3), layout_matrix = lay)
dev.off()

# A&phase&M: power --------------------------------------------
settings = expand.grid(A = c(3, 3.5, 4, 4.5, 5),
                       phase = c(0, 1, 2, 3, 4),
                       M = c(5, 5.5, 6, 6.5, 7),
                       sigma = 1)
settings.keep = sapply(1:nrow(settings), function(i){
  a.settings = settings[i, ]
  rhyth.one.keep =  sum(a.settings$A==3, a.settings$phase ==0, a.settings$M==5, a.settings$sigma==1)>=3
  rhyth.none.keep = a.settings$A==3&a.settings$phase ==0&a.settings$M==5
  return(rhyth.one.keep&(!rhyth.none.keep))
})
settings = settings[settings.keep, ]

time.start = Sys.time()
res = lapply(1:nrow(settings), function(a){
  a.sigma = settings$sigma[a]
  a.M = settings$M[a]
  a.phase = settings$phase[a]
  a.A = settings$A[a]

  if(!file.exists(paste0(out.dir,  "/LR_overall_parameterAPM_simPower", "_A", a.A, "_m", a.M, "_phase", round(a.phase, 2), "_sigma", a.sigma,  ".RData"))){
    B.res = mclapply(1:num_genes, function(b){
      set.seed(b)
      # sample times
      x1.time = runif(n.sample, min=0, max=24)
      x2.time = runif(n.sample, min=0, max=24)

      # noise
      noise.mat1 = matrix(rnorm(1*n.sample, 0, sigma1), ncol=n.sample, nrow=1)
      noise.mat2 = matrix(rnorm(1*n.sample, 0, a.sigma), ncol=n.sample, nrow=1)

      # signal
      signal.mat1 = t(sapply(1:1, function(g){M1 + A1 * cos((2*pi/24)*(x1.time - phase1))}))
      signal.mat2 = t(sapply(1:1, function(g){a.M + a.A * cos((2*pi/24)*(x2.time - a.phase))}))

      x1 = list(data = data.frame(noise.mat1 + signal.mat1),
                time = x1.time,
                gname = paste("gene", seq_len(nrow(signal.mat1))))
      x2 = list(data = data.frame(noise.mat2 + signal.mat2),
                time = x2.time,
                gname = paste("gene", seq_len(nrow(signal.mat2))))
      DCP.obj = DiffCircaPipeline::DCP_Rhythmicity(x1, x2)
      DCP.Dparam = DiffCircaPipeline::DCP_DiffPar(DCP.obj, "A&phase&M", "both")
      testA = diffCircadian::LR_diff(x1.time, noise.mat1 + signal.mat1, x2.time, noise.mat2 + signal.mat2, period = 24, FN = TRUE, type="amplitude")
      testphase = diffCircadian::LR_diff(x1.time, noise.mat1 + signal.mat1, x2.time, noise.mat2 + signal.mat2, period = 24, FN = TRUE, type="phase")
      testM = diffCircadian::LR_diff(x1.time, noise.mat1 + signal.mat1, x2.time, noise.mat2 + signal.mat2, period = 24, FN = TRUE, type="basal")

      a.tab.sidak = data.frame(Combined = DCP.Dparam$p.overall<0.05,
                               A.w.overall = testA$pvalue<0.05, phase.w.overall = testphase$pvalue<0.05,
                               M.w.overall = testM$pvalue<0.05,
                               any.w.P = testA$pvalue<0.05|testphase$pvalue<0.05|testM$pvalue<0.05,
                               A.o.overall = DCP.Dparam$post.hoc.A.By.p,
                               phase.o.overall = DCP.Dparam$post.hoc.peak.By.p,
                               M.o.overall = DCP.Dparam$post.hoc.M.By.p,
                               any.o.overall = DCP.Dparam$post.hoc.A.By.p|DCP.Dparam$post.hoc.peak.By.p|DCP.Dparam$post.hoc.M.By.p,
                               setting = a,
                               method = "Sidak",
                               sigma = settings$sigma[a],
                               M = settings$M[a],
                               phase = settings$phase[a],
                               A = settings$A[a],
                               test.ind = b)

      return(a.tab.sidak)

    }, mc.cores = 40)

    diff.res = do.call(rbind.data.frame, B.res)

    save(diff.res, file=paste0(out.dir,  "/LR_overall_parameterAPM_simPower", "_A", a.A, "_m", a.M, "_phase", round(a.phase, 2), "_sigma", a.sigma,  ".RData"))
  }

  #mark time
  print(paste0(a/nrow(settings)))
  cat("Time spent = "); cat(difftime(Sys.time(), time.start, units = "mins")); cat("\n")
  cat("Time remain = "); cat(difftime(Sys.time(), time.start, units = "mins")/a*(nrow(settings)-a)); cat("\n")

})


# A&phase: power ----------------------------------------------------------
settings = expand.grid(A = c(3, 3.5, 4, 4.5, 5),
                       phase = c(0, 1, 2, 3, 4),
                       M = c(5, 5.5, 6, 6.5, 7),
                       sigma = 1)
settings.keep = sapply(1:nrow(settings), function(i){
  a.settings = settings[i, ]
  rhyth.one.keep =  sum(a.settings$A==3, a.settings$phase ==0, a.settings$M==5, a.settings$sigma==1)>=3
  rhyth.none.keep = a.settings$A==3&a.settings$phase ==0&a.settings$M==5
  return(rhyth.one.keep&(!rhyth.none.keep))
})
settings = settings[settings.keep, ]

time.start = Sys.time()
res = lapply(1:nrow(settings), function(a){
  a.sigma = settings$sigma[a]
  a.M = settings$M[a]
  a.phase = settings$phase[a]
  a.A = settings$A[a]

  if(!file.exists(paste0(out.dir,  "/LR_overall_parameterAP_simPower", "_A", a.A, "_m", a.M, "_phase", round(a.phase, 2), "_sigma", a.sigma,  ".RData"))){
    B.res = mclapply(1:num_genes, function(b){
      set.seed(b)
      # sample times
      x1.time = runif(n.sample, min=0, max=24)
      x2.time = runif(n.sample, min=0, max=24)

      # noise
      noise.mat1 = matrix(rnorm(1*n.sample, 0, sigma1), ncol=n.sample, nrow=1)
      noise.mat2 = matrix(rnorm(1*n.sample, 0, a.sigma), ncol=n.sample, nrow=1)

      # signal
      signal.mat1 = t(sapply(1:1, function(g){M1 + A1 * cos((2*pi/24)*(x1.time - phase1))}))
      signal.mat2 = t(sapply(1:1, function(g){a.M + a.A * cos((2*pi/24)*(x2.time - a.phase))}))

      x1 = list(data = data.frame(noise.mat1 + signal.mat1),
                time = x1.time,
                gname = paste("gene", seq_len(nrow(signal.mat1))))
      x2 = list(data = data.frame(noise.mat2 + signal.mat2),
                time = x2.time,
                gname = paste("gene", seq_len(nrow(signal.mat2))))
      DCP.obj = DiffCircaPipeline::DCP_Rhythmicity(x1, x2)
      DCP.Dparam = DiffCircaPipeline::DCP_DiffPar(DCP.obj, "A&phase", "both")
      testA = diffCircadian::LR_diff(x1.time, noise.mat1 + signal.mat1, x2.time, noise.mat2 + signal.mat2, period = 24, FN = TRUE, type="amplitude")
      testphase = diffCircadian::LR_diff(x1.time, noise.mat1 + signal.mat1, x2.time, noise.mat2 + signal.mat2, period = 24, FN = TRUE, type="phase")
      # testM = diffCircadian::LR_diff(x1.time, noise.mat1 + signal.mat1, x2.time, noise.mat2 + signal.mat2, period = 24, FN = TRUE, type="basal")

      a.tab.sidak = data.frame(Combined = DCP.Dparam$p.overall<0.05,
                               A.w.overall = testA$pvalue<0.05, phase.w.overall = testphase$pvalue<0.05,
                               # M.w.overall = testM$pvalue<0.05,
                               any.w.P = testA$pvalue<0.05|testphase$pvalue<0.05,
                               A.o.overall = DCP.Dparam$post.hoc.A.By.p,
                               phase.o.overall = DCP.Dparam$post.hoc.peak.By.p,
                               # M.o.overall = DCP.Dparam$post.hoc.M.By.p,
                               any.o.overall = DCP.Dparam$post.hoc.A.By.p|DCP.Dparam$post.hoc.peak.By.p,
                               setting = a,
                               method = "Sidak",
                               sigma = settings$sigma[a],
                               M = settings$M[a],
                               phase = settings$phase[a],
                               A = settings$A[a],
                               test.ind = b)

      return(a.tab.sidak)

    }, mc.cores = 40)

    diff.res = do.call(rbind.data.frame, B.res)

    save(diff.res, file=paste0(out.dir,  "/LR_overall_parameterAP_simPower", "_A", a.A, "_m", a.M, "_phase", round(a.phase, 2), "_sigma", a.sigma,  ".RData"))
  }

  #mark time
  print(paste0(a/nrow(settings)))
  cat("Time spent = "); cat(difftime(Sys.time(), time.start, units = "mins")); cat("\n")
  cat("Time remain = "); cat(difftime(Sys.time(), time.start, units = "mins")/a*(nrow(settings)-a)); cat("\n")

})

# Summary -----------------------------------------------------------------
res = mclapply(1:nrow(settings), function(a){
  a.sigma = settings$sigma[a]
  a.M = settings$M[a]
  a.phase = settings$phase[a]
  a.A = settings$A[a]

  load(file=paste0(out.dir,  "/LR_overall_parameterAPM_simPower", "_A", a.A, "_m", a.M, "_phase", round(a.phase, 2), "_sigma", a.sigma,  ".RData"))
  return(diff.res)
  #mark time
},mc.cores = 10)

rep = rep(1:10, each = 10000)
res = lapply(res, function(a.res){
  a.res$rep = rep
  return(a.res)
})

res.tab = do.call(rbind.data.frame, res)
res.tab.mean.s = res.tab %>%
  split(f = list(res.tab$M, res.tab$phase, res.tab$A, res.tab$sigma, res.tab$rep, res.tab$method)) %>%
  lapply(function(x){
    a.tab = data.frame(Combined = mean(x$Combined),
                       # A.alone = mean(x$A.w.overall), phase.alone = mean(x$phase.w.overall), M.alone = mean(x$M.w.overall),
                       A.o.overall = mean(x$A.o.overall), phase.o.overall = mean(x$phase.o.overall), M.o.overall = mean(x$M.o.overall),
                       any.o.overall = mean(x$any.o.overall),
                       method = x$method[1],
                       sigma = x$sigma[1],
                       M2 = x$M[1],
                       phase2 = x$phase[1],
                       A2 = x$A[1])
    return(a.tab)
  })
res.tab.mean = do.call(rbind.data.frame, res.tab.mean.s)
res.tab.mean = res.tab.mean[complete.cases(res.tab.mean), ]

res.tab.mean = res.tab.mean %>% pivot_longer(cols = c(1:5),
                                             names_to = "Tests",
                                             values_to = "Power")
res.tab.mean$Tests = sapply(res.tab.mean$Tests, function(a){
  switch(a,
         "Combined" = "Combined",
         "A.o.overall" = "Post hoc A",
         "phase.o.overall" = "Post hoc phase",
         "M.o.overall" = "Post hoc M",
         "any.o.overall" = "Post hoc any")
})
res.tab.mean$Tests = factor(res.tab.mean$Tests, levels = c("Combined", "Post hoc A", "Post hoc phase", "Post hoc M", "Post hoc any"))

#plot with change of A
res.to.A = res.tab.mean[res.tab.mean$M2==5&res.tab.mean$phase2==0&res.tab.mean$sigma==1, ]
res.to.A = res.to.A %>%
  split(f = list(res.to.A$A2, res.to.A$Tests)) %>%
  lapply(function(x){
    if(nrow(x)==0){
      return(NULL)
    }else{
      x.tab = data.frame(sigma = x$sigma[1],
                         M = x$M2[1],
                         phase = x$phase2[1],
                         A = x$A2[1],
                         mean = mean(x$Power),
                         sd = sd(x$Power),
                         Tests = x$Tests[1])
      return(x.tab)
    }
  })
res.to.A = do.call(rbind.data.frame, res.to.A)
p3.pwer.1 = ggplot(res.to.A, aes(x = A, y = mean, color = Tests))+
  geom_line(position = position_dodge(0.2))+
  geom_point(alpha = 0.8, size = 3, position = position_dodge(0.2))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.4,
                position=position_dodge(.2))+
  # scale_y_continuous(breaks=seq(0, 0.35, by = 0.02))+
  theme_bw()+
  scale_color_manual(labels = c(expression(Overall~(Delta*A*","~Delta*phi*","~Delta*M)), expression(Post*"-"*hoc~Delta*A),
                               expression(Post*"-"*hoc~Delta*phi), expression(Post*"-"*hoc~Delta*M),
                               expression(Post*"-"*hoc~(Delta*A*","~Delta*phi*","~Delta*M))),
                    values = c("#F8766D", "#D39200", "#00BA38", "#00B9E3", "#DB72FB"))+
  labs(title = expression(C*.~Power~with~change~of~A[2]~(A[1]~"="~3)), x = expression(A[2]), y = "Power")+
  theme(legend.position = "none", legend.text.align = 0)

res.to.phase = res.tab.mean[res.tab.mean$A2==3&res.tab.mean$M2==5&res.tab.mean$sigma==1&res.tab.mean$phase2!=0, ]
res.to.phase = res.to.phase %>%
  split(f = list(res.to.phase$phase2, res.to.phase$Tests)) %>%
  lapply(function(x){
    if(nrow(x)==0){
      return(NULL)
    }else{
      x.tab = data.frame(sigma = x$sigma[1],
                         M = x$M2[1],
                         phase = x$phase2[1],
                         A = x$A2[1],
                         mean = mean(x$Power),
                         sd = sd(x$Power),
                         Tests = x$Tests[1])
      return(x.tab)
    }
  })
res.to.phase = do.call(rbind.data.frame, res.to.phase)
p3.pwer.2 = ggplot(res.to.phase, aes(x = phase, y = mean, color = Tests))+
  geom_line(position = position_dodge(0.2))+
  geom_point(alpha = 0.8, size = 3, position = position_dodge(0.2))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.2))+
  scale_x_continuous(breaks=c(1, 2, 3, 4), labels = c(1, 2, 3, 4))+
  scale_color_manual(labels = c(expression(Overall~(Delta*A*","~Delta*phi*","~Delta*M)),
                                expression(Post*"-"*hoc~Delta*A),
                                expression(Post*"-"*hoc~Delta*phi),
                                expression(Post*"-"*hoc~Delta*M),
                                expression(Post*"-"*hoc~(Delta*A*","~Delta*phi*","~Delta*M))),
                     values = c("#F8766D", "#D39200", "#00BA38", "#00B9E3", "#DB72FB"))+
  # scale_x_continuous(breaks=seq(0, 0.35, by = 0.02)s)+0, pi/16, pi/8,pi/16*3, pi/4
  theme_bw()+
  labs(title = expression(D*.~Power~with~change~of~phi[2]~(phi[1]~"="~0)), x = expression(phi[2]), y = "Power")+
  theme(legend.position = "none", legend.text.align = 0)

res.to.M = res.tab.mean[res.tab.mean$A2==3&res.tab.mean$phase2==0&res.tab.mean$sigma==1&res.tab.mean$M2!=5, ]
res.to.M = res.to.M %>%
  split(f = list(res.to.M$M2, res.to.M$Tests)) %>%
  lapply(function(x){
    if(nrow(x)==0){
      return(NULL)
    }else{
      x.tab = data.frame(sigma = x$sigma[1],
                         M = x$M2[1],
                         phase = x$phase2[1],
                         A = x$A2[1],
                         mean = mean(x$Power),
                         sd = sd(x$Power),
                         Tests = x$Tests[1])
      return(x.tab)
    }
  })
res.to.M = do.call(rbind.data.frame, res.to.M)
p3.pwer.3 = ggplot(res.to.M, aes(x = M, y = mean, color = Tests))+
  geom_line(position = position_dodge(0.2))+
  geom_point(alpha = 0.8, size = 3, position = position_dodge(0.2))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.4,
                position=position_dodge(.2))+
  scale_color_manual(labels = c(expression(Overall~(Delta*A*","~Delta*phi*","~Delta*M)), expression(Post*"-"*hoc~Delta*A),
                                expression(Post*"-"*hoc~Delta*phi), expression(Post*"-"*hoc~Delta*M),
                                expression(Post*"-"*hoc~(Delta*A*","~Delta*phi*","~Delta*M))),
                     values = c("#F8766D", "#D39200", "#00BA38", "#00B9E3", "#DB72FB"))+
  theme_bw()+
  labs(title = expression(E*.~Power~with~change~of~M[2]~(M[1]~"="~5)), x = expression(M[2]), y = "Power")+
  theme(legend.text.align = 0)

# pdf(paste0(out.dir, "/Paper_LR_overall_parameter3_Power_all.pdf"), width = 15, height = 5)
# lay <- matrix(c(1,1,1,2,2,2, 3, 3,3,3), ncol = 10)
# gridExtra::grid.arrange(grobs = list(p3.pwer.1, p3.pwer.2, p3.pwer.3), layout_matrix = lay)
# dev.off()

#plot power for A&phase
settings = expand.grid(A = c(3, 3.5, 4, 4.5, 5),
                       phase = c(0, 1, 2, 3, 4),
                       M = c(5, 5.5, 6, 6.5, 7),
                       sigma = 1)
settings.keep = sapply(1:nrow(settings), function(i){
  a.settings = settings[i, ]
  rhyth.one.keep =  sum(a.settings$A==3, a.settings$phase ==0, a.settings$M==5, a.settings$sigma==1)>=3
  rhyth.none.keep = a.settings$A==3&a.settings$phase ==0&a.settings$M==5
  return(rhyth.one.keep&(!rhyth.none.keep))
})
settings = settings[settings.keep, ]

res = mclapply(1:nrow(settings), function(a){
  a.sigma = settings$sigma[a]
  a.M = settings$M[a]
  a.phase = settings$phase[a]
  a.A = settings$A[a]

  load(file=paste0(out.dir,  "/LR_overall_parameterAP_simPower", "_A", a.A, "_m", a.M, "_phase", round(a.phase, 2), "_sigma", a.sigma,  ".RData"))
  return(diff.res)
  #mark time
},mc.cores = 10)

rep = rep(1:10, each = 10000)
res = lapply(res, function(a.res){
  a.res$rep = rep
  return(a.res)
})

res.tab = do.call(rbind.data.frame, res)
res.tab.mean.s = res.tab %>%
  split(f = list(res.tab$setting, res.tab$rep, res.tab$method)) %>%
  lapply(function(x){
    a.tab = data.frame(Combined = mean(x$Combined),
                       A.o.overall = mean(x$A.o.overall), phase.o.overall = mean(x$phase.o.overall), any.o.overall = mean(x$any.o.overall),
                       method = x$method[1],
                       sigma = x$sigma[1],
                       M = x$M[1],
                       phase = x$phase[1],
                       A = x$A[1])
    return(a.tab)
  })
res.tab.mean = do.call(rbind.data.frame, res.tab.mean.s)
# res.tab.mean = res.tab.mean[complete.cases(res.tab.mean), ]
res.tab.mean$M = factor(res.tab.mean$M)
res.tab.mean$A = factor(res.tab.mean$A)
res.tab.mean$phase = factor(round(res.tab.mean$phase,2))
res.tab.mean$sigma = factor(round(res.tab.mean$sigma,2))

res.tab.mean = res.tab.mean %>% pivot_longer(cols = c(1:4),
                                             names_to = "Tests",
                                             values_to = "Power")
res.tab.mean$Tests = sapply(res.tab.mean$Tests, function(a){
  switch(a,
         "any.o.overall" = "Post hoc any",
         "A.o.overall" = "Post hoc A",
         "phase.o.overall" = "Post hoc phase",
         "Combined" = "Combined")
})
res.tab.mean$Tests = factor(res.tab.mean$Tests, levels = c("Combined", "Post hoc A", "Post hoc phase", "Post hoc any"))

res.to.A = res.tab.mean[res.tab.mean$M==5&res.tab.mean$phase==0&res.tab.mean$sigma==1&res.tab.mean$A!=3, ]
res.to.A = res.to.A %>%
  split(f = list(res.to.A$A, res.to.A$Tests)) %>%
  lapply(function(x){
    if(nrow(x)==0){
      return(NULL)
    }else{
      x.tab = data.frame(sigma = x$sigma[1],
                         M = x$M[1],
                         phase = x$phase[1],
                         A = x$A[1],
                         mean = mean(x$Power),
                         sd = sd(x$Power),
                         Tests = x$Tests[1])
      return(x.tab)
    }
  })
res.to.A = do.call(rbind.data.frame, res.to.A)

res.to.phase = res.tab.mean[res.tab.mean$A==3&res.tab.mean$M==5&res.tab.mean$sigma==1&res.tab.mean$phase!=0, ]
res.to.phase = res.to.phase %>%
  split(f = list(res.to.phase$phase, res.to.phase$Tests)) %>%
  lapply(function(x){
    if(nrow(x)==0){
      return(NULL)
    }else{
      x.tab = data.frame(sigma = x$sigma[1],
                         M = x$M[1],
                         phase = x$phase[1],
                         A = x$A[1],
                         mean = mean(x$Power),
                         sd = sd(x$Power),
                         Tests = x$Tests[1])
      return(x.tab)
    }
  })
res.to.phase = do.call(rbind.data.frame, res.to.phase)

p2.pwer.1 = ggplot(res.to.A %>% mutate(A = as.numeric(A)), aes(x = A, y = mean, color = Tests))+
  geom_line(position = position_dodge(0.2))+
  geom_point(alpha = 0.8, size = 3, position = position_dodge(0.2))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.4,
                position=position_dodge(.2))+
  theme_bw()+
  scale_color_manual(labels = c(expression(Combined~(Delta*A*","~Delta*phi)), expression(Post*"-"*hoc~Delta*A),
                                expression(Post*"-"*hoc~Delta*phi),
                                expression(Post*"-"*hoc~(Delta*A*","~Delta*phi))),
                     values = c("#F8766D", "#D39200", "#00BA38", "#00B9E3"))+
  labs(title = expression(A*.~Power~with~change~of~A[2]~(A[1]~"="~3)), x = expression(A[2]), y = "Power")+
  theme(legend.position = "none", legend.text.align = 0)

p2.pwer.2 = ggplot(res.to.phase %>% mutate(phase = as.numeric(as.character(phase))), aes(x = phase, y = mean, color = Tests))+
  geom_line(position = position_dodge(0.2))+
  geom_point(alpha = 0.8, size = 3, position = position_dodge(0.2))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.2))+
  scale_x_continuous(breaks=c(1, 2, 3, 4), labels = c(1, 2, 3, 4))+
  scale_color_manual(labels = c(expression(Combined~(Delta*A*","~Delta*phi)), expression(Post*"-"*hoc~Delta*A),
                                expression(Post*"-"*hoc~Delta*phi),
                                expression(Post*"-"*hoc~(Delta*A*","~Delta*phi))),
                     values = c("#F8766D", "#D39200", "#00BA38", "#00B9E3"))+
  # scale_x_continuous(breaks=seq(0, 0.35, by = 0.02)s)+0, pi/16, pi/8,pi/16*3, pi/4
  theme_bw()+
  labs(title = expression(B*.~Power~with~change~of~phi[2]~(phi[1]~"="~0)), x = expression(phi[2]), y = "Power")+
  theme( legend.text.align = 0)

pdf(paste0(out.dir, "/paper_LR_overall_Power.pdf"), width = 15, height = 10)
lay <- matrix(
  c(6,1,1,1,2,2,2,2,6,6,
    3,3,3,4,4,4,5,5,5,5), ncol = 10, byrow = TRUE)
gridExtra::grid.arrange(grobs = list(p2.pwer.1, p2.pwer.2,
                                     p3.pwer.1, p3.pwer.2, p3.pwer.3), layout_matrix = lay)
dev.off()















