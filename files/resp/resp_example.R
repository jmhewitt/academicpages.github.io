## ----install, cache=F----------------------------------------------------

# simple looping
if (!require("foreach")) install.packages("foreach")

# parallelization
if (!require("doMC")) install.packages("doMC")
if (!require("doRNG")) install.packages("doRNG")

# data manipulation and plotting
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("cowplot")) install.packages("cowplot")

# MCMC tools
if (!require("coda")) install.packages("coda")

# RESP model
if (!require("devtools")) install.packages("devtools")
if (!require("telefit")) devtools::install_github('jmhewitt/telefit')

## ----load----------------------------------------------------------------
load('resp_coex.RData')
str(dat)

## ----fit, eval=F---------------------------------------------------------
## 
## # sampler parameters (Note: cov.r "nugget" is not currently estimated)
## maxIt = 21000
## burn = 1000
## priors = list(
##   beta = list( Lambda = diag(10, ncol(dat$X[,,1]))),
##   cov.s = list( smoothness = .5, range = c(1, 600), variance = c(2, 30),
##                 nugget = c(2, 1) ),
##   cov.r = list( smoothness = .5, range = c(1, 600), variance = c(2,1e-1),
##                 nugget = c(2,1) )
## )
## 
## # parallelize estimation of (sub)models
## ncores = detectCores() - 1
## registerDoMC(ncores)
## 
## # fit (sub)models in parallel
## ret = foreach(i=1:3) %dorng% {
##   if(i==1) { # RESP model
##     stFit(dat, priors, maxIt, coords.knots = coords.knots)
##   } else if (i==2) { # SP submodel
##     stFit(dat, priors, maxIt, coords.knots = coords.knots, localOnly = T)
##   } else if (i==3) { # RE submodel
##     stFit(dat, priors, maxIt, coords.knots = coords.knots, remoteOnly = T)
##   }
## }
## 
## # extract model fits
## fit = ret[[1]]
## fit.local = ret[[2]]
## fit.remote = ret[[3]]
## rm(ret)

## ----posterior_samples---------------------------------------------------

# posterior parameter samples
str(fit$parameters$samples)

## ----composition_sampling, eval=F----------------------------------------
## 
## # estimate teleconnection effects and make predictions
## fcst.local = stPredict(stFit = fit.local, stData = dat,
##                        stDataNew = dat, burn = burn, ncores = ncores)
## fcst.remote = stPredict(stFit = fit.remote, stData = dat, stDataNew = dat,
##                         burn = burn, ncores = ncores, returnFullAlphas = F)
## fcst = stPredict(stFit = fit, stData = dat, stDataNew = dat, burn = burn,
##                  ncores = ncores, returnFullAlphas = T)

## ----evaluate_predictions, results='markup'------------------------------

# evaluate predictions
clim = rowMeans(dat$Y)
fcst = stEval(fcst, dat$Y, clim)
fcst.local = stEval(fcst.local, dat$Y, clim)
fcst.remote = stEval(fcst.remote, dat$Y, clim)

# look at errors for predictions at first timepoint
str(fcst$pred[[1]]$err)

# compute VIFs
vif = stVIF(stData = dat, stFit = fit, burn = 1000)
vif$beta
summary(vif$alpha)

## ----figure_1, fig.height=3, fig.width=8---------------------------------

# unnormalize data
dat$Z = dat$Z*5252

# set plotting label
dat$Z.lab = 'z'

# set plot year
t = 1982

# build base plot of sea surface temperatures
x = plot(dat, type='remote', t=t, boxsize = 6.5) + ggtitle('')

# re-normalize data
dat$Z = dat$Z/5252

# add local precipitation
x = x + 
  geom_point(aes(x=lon, y=lat, colour = precip), inherit.aes = F, size = .6,
             stroke = 0, pch = 15,
             data = data.frame(lon = dat$coords.s[,1], lat = dat$coords.s[,2],
                               precip = dat$Y[,which(dat$tLabs==t)])) +
  scale_color_distiller('Y', palette = 'PuBuGn', direction = 1) +
  theme(legend.box = 'horizontal',
        legend.box.margin = margin(t=0, r=-8, b=0, l=0, unit='pt'))

# set left endpoints of arrows for schematic
arrows.start = data.frame(
  x = c(-100, -150, -200, -180),
  y = c(-5, 5, 35, 55)
)[c(3,2,1,4),] %>% mutate(x=ifelse(x<=0,x,x-360))

# set arrow color
col = 'forestgreen'

# plot after adding schematic arrows and labels
x + geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),
               curvature = .3, arrow = arrow(length = unit(0.03, "npc")),
               inherit.aes = F, lwd = 1, col=col,
               data = data.frame(x1=arrows.start$x[1:3],
                                 y1=arrows.start$y[1:3],
                                 x2=dat$coords.s[c(161,187,219),1],
                                 y2=dat$coords.s[c(161,187,219),2])
               ) +
  geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),
             curvature = -.3, arrow = arrow(length = unit(0.03, "npc")),
             inherit.aes = F, lwd = 1, col=col,
             data = data.frame(x1=arrows.start$x[4],
                               y1=arrows.start$y[4],
                               x2=dat$coords.s[23,1],
                               y2=dat$coords.s[23,2])
             ) +
  geom_point(aes(x=x, y=y), size=3, col=col,
             data = arrows.start, inherit.aes = F, pch=19) +
  geom_label(aes(x=x,y=y,label=label), inherit.aes = F,
             data = data.frame(
               x=c(-185,-110),
               y=c(15,dat$coords.s[94,2]+10),
               label=c('z(r,t)','x(s,t), Y(s,t)')
             ), col=col, size=5)

## ----figure_2, fig.width=8, fig.height=3---------------------------------

plot_grid(
  plot(dat, type='eof') +
    ggtitle('SST Empirical orthogonal function 1 (EOF 1)') +
    theme(text = element_text(size=8)),
  plot(dat, type='eof_cor', signif.level = .05, signif.telecon = T,
       lwd=.7, alpha=.6) +
    ggtitle('Correlation between Precip. and EOF 1 score') +
    theme(text = element_text(size=8)),
  ncol=2,
  labels=paste(LETTERS[1:2], ')', sep='')
)

## ----figure_3, fig.width=8, fig.height=3---------------------------------

# build base plot of teleconnection effect estimates
resp = plot(fcst, type='eof_alpha_knots', stFit = fit, stData = dat,
            signif.telecon = T, lwd=0.7, alpha=.6) +
  scale_fill_gradient2(expression(hat(alpha)~"'"(~bold(".")~","~1)),
                       low = "#0571b0", mid = '#f7f7f7', high = '#ca0020') +
  ggtitle('Estimated teleconnection effects associated with EOF 1') 

# plot after adding denver as a point of reference
resp +
  geom_point(aes(x=lon, y=lat), 
             data = data.frame(lon = -104.9903, lat = 39.7392),
             inherit.aes = F, size = 2) + 
  geom_label(aes(x=lon, y=lat, label=label),
             data = data.frame(lon = -103.5, lat = 39.7392, label = 'Denver'),
             inherit.aes = F, alpha = .85) + 
  theme(text = element_text(size=8))

## ----figure_4, fig.width=7, fig.height=3---------------------------------

# set grouping labels
ilabs = c('Submodels', 'Common models')

# build plot data
errs = errs %>% mutate(model = factor(model, levels = c('RESP', 'RE', 'SP', 
                                                        'SVC', 'ENSO-T', 'CCA', 
                                                        'CLIM')),
         crps.rel = crps.cat / median(errs$crps.cat[errs$model=='CLIM']),
         comp = recode(model, 'RE'=ilabs[1], 'SP'=ilabs[1], 'RESP'='',
                       'SVC'=ilabs[2], 'ENSO-T'=ilabs[2], 'CCA'=ilabs[2], 
                       'CLIM'=ilabs[2]))

# plot error comparison
errs %>% ggplot(aes(x=model, y=crps.rel)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, lty = 2) +
  xlab('') +
  ylab('Relative RPS') +
  scale_y_continuous(breaks = c(.5, 1, 1.5, 2, 3)) +
  facet_grid(~comp, scales = 'free_x', space = 'free_x') +
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = 'gray15'),
        panel.border = element_rect(colour = 'gray60', linetype = 1, size = .5),
        axis.line.y = element_blank(),
        panel.spacing = unit(0, 'pt'))

## ----supplement_figure_c1, fig.width=7, fig.height=3---------------------

errs %>% ggplot(aes(x=model, y=cat.heidke.alt)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty=2) +
  xlab('') +
  ylab('Heidke skill score') +
  coord_equal(ratio = 2) +
  facet_grid(~comp, scales = 'free_x', space = 'free_x') +
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = 'gray15'),
        panel.border = element_rect(colour = 'gray60', linetype = 1, size = .5),
        axis.line.y = element_blank(),
        panel.spacing = unit(0, 'pt'))

## ----table_1-------------------------------------------------------------

maxIt = 21e3
burn = 1e3

# extract posterior samples of parameters
params = data.frame(fit$parameters$samples) %>% 
  select(-ll, -sigmasq_r_eps) %>%
  slice(burn:maxIt)
colnames(params)[1:ncol(dat$X)] = colnames(dat$X)

# "unexpand" the nugget parameter
params[,5] = params[,5] * params[,3]

# compute highest posterior density intervals
hpds = HPDinterval(mcmc(params))
hpds = round(hpds, 2)

# compute posterior estimates
ests = round(colMeans(mcmc(params)),2)

# stronger rounding for range parameters
hpds[6:7,] = round(hpds[6:7,], 0)
ests[6:7] = round(ests[6:7], 0)

# assemble parameter estimates table
df = data.frame(
  Param = c('$\\beta_0$', '$\\beta_T$', '$\\sigma^2_w$', '$\\sigma^2_\\alpha$',
            '$\\sigma^2_\\varepsilon$', '$\\rho_w$', '$\\rho_\\alpha$'),
  Est=ests,
  HPD=paste('(', paste(hpds[,1], hpds[,2], sep=', '), ')', sep='')
)

# text formatting
df$Est = gsub('-', '$-$',sprintf('%.2f', df$Est))
df[,1] = as.character(df[,1])
df[,3] = gsub('-','$-$',as.character(df[,3]))

colnames(df) = c('', 'Posterior mean', '95\\% HPD')
rownames(df) = c('Local effects',' ',  '  ', '   ',
                 'Covariance', '    ', '     ')

df

## ----supplement_figure_d1, fig.width=8, fig.height=6---------------------

plot_grid(
  plot(dat.test, fill.lab=expression(Y)) + ggtitle('Observed (1982)'),
  plot(fcst, fill.lab=expression(hat(Y))) + ggtitle('Forecast') +
    theme(axis.title.y = element_text(color='white')),
  NULL,
  plot(fcst, type='se', fill.lab=expression(SE(hat(Y)))) + 
    ggtitle('Uncertainty'),
  ncol=2, labels = c('A)','B)', '', 'C)')
)

## ----supplement_figure_d2, fig.width=8, fig.height=6---------------------

# compute posterior logits
logits = qlogis(fcst.test$samples$cat_probs)
logits.est = matrix(nrow = nrow(logits), ncol = ncol(logits))
for(i in 1:nrow(logits)) {
  logits.est[i,] = rowMeans(logits[i,,])
}

# transfer logits to a plottable object
logitDat = dat.test
logitDat$Y = apply(logits.est, 1, max) 

# build plots

p = plot(dat.test, type='cat.response', 
         category.breaks = fcst$category.breaks, fill.lab='') +
  ggtitle('Observed (1982)')

uncty = plot(logitDat, fill.lab=expression(hat(logit)(tilde(Y)))) +
  ggtitle('Uncertainty')

# assemble and display plots
ggdraw(plot_grid(
  plot_grid(
    plot_grid(
      p + theme(legend.position = 'none'),
      NULL,
      plot(fcst.test, type='cat.pred', fill.lab='') +
        ggtitle('Forecast') +
        theme(axis.title.y = element_text(color='white'),
              legend.position = 'none'),
      nrow=1, align='v', rel_widths = c(1,.03,1), labels=c('A)','','B)')),
    get_legend(p),
  rel_widths=c(1, 0.2)
  ),
    plot_grid(
    plot_grid(
      NULL,
      NULL,
      uncty + theme(legend.position = 'none'),
      nrow=1, align='v', rel_widths = c(1,.03,1), labels=c('','','C)')),
    get_legend(uncty),
  rel_widths=c(1, 0.2)
  ),
  nrow=2
))

## ----supplement_figure_d3, fig.width=8, fig.height=3---------------------

dat$Z.lab = 'SST'
plot(dat, type='remote', coords.knots = coords.knots, t = 1982) +
  ggtitle('')

