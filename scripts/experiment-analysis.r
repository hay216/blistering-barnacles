


library(lme4)    # For mixed-effects models

library(equivalence) # For data that I know

library(ggplot2)

library(reshape2)

data(ufc)

str(ufc)

ufc <- subset(ufc, !is.na(Dbh))
ufc$Species[ufc$Species %in% c("F","FG")] <- "GF"
ufc$Species <- factor(ufc$Species)
ufc$ht.measured <- is.na(ufc$Height)

ufc <- subset(ufc,
              Species %in% names(sort(-table(Species)))[1:6],
              select = c("ht.measured","Plot","Dbh.in","Species"))

## Fit a pair of models to compare using LRT

ufc.whole <- ufc

isna.1.w <- glmer(ht.measured ~ Dbh.in + Species + (1 | Plot),
                  data = ufc.whole, family = binomial)
isna.0.w <- glmer(ht.measured ~ Dbh.in + (1 | Plot),
                  data = ufc.whole, family = binomial)

anova(isna.0.w, isna.1.w)

load("../images/mi.RData")

ls()

str(experiment)

experiment$naive.p <- unlist(experiment$naive.p)
experiment$imputed.p <- unlist(experiment$imputed.p)

experiment$ms <- factor(experiment$ms)
experiment$ns <- factor(experiment$ns)

levels(experiment$ms) <- paste("m = ", levels(experiment$ms))
levels(experiment$ns) <- paste("n = ", levels(experiment$ns))

qplot(x = naive.p,
      geom = "density",
      data = experiment) + facet_grid(ms ~ ns, scales = "free") +
  theme(axis.text.x = element_text(angle = 45)) 

qplot(x = imputed.p,
      geom = "density",
      data = experiment) + facet_grid(ms ~ ns) +
  theme(axis.text.x = element_text(angle = 45)) 

exp.melt <- melt(experiment,
                 id.vars = c("ms","ns","reps"),
                 variable.name = "Estimate")

head(exp.melt)

qplot(x = value,
      geom = "density",
      fill = variable,
      data = exp.melt) + facet_grid(ms ~ ns) +
  theme(axis.text.x = element_text(angle = 45)) 

