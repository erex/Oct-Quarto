#| label = "setup",
#| include = FALSE
knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning = FALSE)
options(width=90)



library(Distance)
data(sikadeer)
conversion.factor <- convert_units("centimeter", "kilometer", "square kilometer")



deer.df <- ds(sikadeer, key="hn", truncation="10%", convert_units = conversion.factor)
plot(deer.df)
print(deer.df$dht$individuals$summary)


#| echo = FALSE,
#| fig.height = 4
# Calculate dung decay rate parameters
MIKE.persistence <- function(DATA) {
  
#  Purpose: calculate mean persistence time (mean time to decay) for dung/nest data 
#  Input: data frame with at least two columns:
#         DAYS - calendar day on which dung status was observed
#         STATE - dung status: 1-intact, 0-decayed
#  Output: point estimate, standard error and CV of mean persistence time
#
#  Attribution: code from Mike Meredith website: 
#      http://www.mikemeredith.net/blog/2017/Sign_persistence.htm
#   Citing: CITES elephant protocol
#      https://cites.org/sites/default/files/common/prog/mike/survey/dung_standards.pdf
  
  ##   Fit logistic regression model to STATE on DAYS, extract coefficients
  dung.glm <- glm(STATE ~ DAYS, data=DATA, family=binomial(link = "logit"))
  betas <- coefficients(dung.glm)
  
  ##   Calculate mean persistence time
  mean.decay <- -(1+exp(-betas[1])) * log(1+exp(betas[1])) / betas[2]
  
  ## Calculate the variance of the estimate
  vcovar <- vcov(dung.glm)
  var0 <- vcovar[1,1]  # variance of beta0
  var1 <- vcovar[2,2]  # variance of beta1
  covar <- vcovar[2,1] # covariance
  
  deriv0 <- -(1-exp(-betas[1]) * log(1+exp(betas[1])))/betas[2]
  deriv1 <- -mean.decay/betas[2]
  
  var.mean <- var0*deriv0^2 + 2*covar*deriv0*deriv1 + var1*deriv1^2
  
  ## Calculate the SE and CV and return
  se.mean <- sqrt(var.mean)
  cv.mean <- se.mean/mean.decay
  
  out <- c(mean.decay, se.mean, 100*cv.mean)
  names(out) <- c("Mean persistence time", "SE", "%CV")
  plot(decay$DAYS, jitter(decay$STATE, amount=0.10), xlab="Days since initiation",
       ylab="Dung persists (yes=1)",
       main="Eight dung piles revisited over time")
  curve(predict(dung.glm, data.frame(DAYS=x), type="resp"), add=TRUE)
  abline(v=mean.decay, lwd=2, lty=3)
  return(out)
}
decay <- read.csv(file="https://raw.githubusercontent.com/erex/Oct-Quarto/main/Pr9/IntroDS_9.1.csv")
persistence.time <- MIKE.persistence(decay)
print(persistence.time)



# Create list of multipliers
mult <- list(creation = data.frame(rate=25, SE=0),
             decay    = data.frame(rate=163, SE=14))
print(mult)
deer_ests <- dht2(deer.df, flatfile=sikadeer, strat_formula=~Region.Label,
                  convert_units=conversion.factor, multipliers=mult, 
                  stratification="effort_sum", total_area = 100)
print(deer_ests, report="abundance")



data(CueCountingExample)
head(CueCountingExample, n=3)
CueCountingExample$Effort <- CueCountingExample$Search.time
cuerates <- CueCountingExample[ ,c("Cue.rate", "Cue.rate.SE", "Cue.rate.df")]
cuerates <- unique(cuerates)
names(cuerates) <- c("rate", "SE", "df")
mult <- list(creation=cuerates)
print(mult)



# Tidy up data by removing surplus columns
CueCountingExample[ ,c("Cue.rate", "Cue.rate.SE", "Cue.rate.df", "Sample.Fraction", 
                       "Sample.Fraction.SE")] <- list(NULL)
trunc <- 1.2
whale.df.hn <- ds(CueCountingExample, key="hn", transect="point", adjustment=NULL,
                  truncation=trunc)
whale.df.hr <- ds(CueCountingExample, key="hr", transect="point", adjustment=NULL,
                  truncation=trunc)
knitr::kable(summarize_ds_models(whale.df.hn, whale.df.hr), digits = 3)



whale.est.hn <- dht2(whale.df.hn, flatfile=CueCountingExample, strat_formula=~1, 
                     multipliers=mult, sample_fraction=0.5)
print(whale.est.hn, report="abundance")
whale.est.hr <- dht2(whale.df.hr, flatfile=CueCountingExample, strat_formula=~1, 
                     multipliers=mult, sample_fraction=0.5)
print(whale.est.hr, report="abundance")


#| layout-ncol: 2
#| fig.cap:
#|   - "Half normal PDF"
#|   - "Hazard rate PDF"
plot(whale.df.hn, pdf=TRUE, main="Half normal")
plot(whale.df.hr, pdf=TRUE, main="Hazard rate")



data(wren_cuecount)
cuerate <- unique(wren_cuecount[ , c("Cue.rate","Cue.rate.SE")])
names(cuerate) <- c("rate", "SE")
mult <- list(creation=cuerate)
print(mult)
# Search time is the effort - this is 2 * 5min visits
wren_cuecount$Effort <- wren_cuecount$Search.time
w3.hr <- ds(wren_cuecount, transect="point", key="hr", adjustment=NULL, truncation=92.5)



conversion.factor <- convert_units("meter", NULL, "hectare")
w3.est <- dht2(w3.hr, flatfile=wren_cuecount, strat_formula=~1,
               multipliers=mult, convert_units=conversion.factor)
# NB "Effort" here is sum(Search.time) in minutes
# NB "CoveredArea" here is pi * w^2 * sum(Search.time)
print(w3.est, report="density")


#| layout-ncol: 2
plot(w3.hr, pdf=TRUE, main="Cue distances of winter wren.")
gof_ds(w3.hr)

