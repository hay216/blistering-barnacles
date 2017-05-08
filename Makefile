# Time-stamp: <2017-05-08 12:39:15 (slane)>
.PHONY: all models input-data output-data clean-models clean-manuscripts clobber \
	PROC-DATA processed-data

all: manuscripts/censored-mle.html \
	processed-data \
	manuscripts/model-interrogation.html \
	manuscripts/censored-mle.pdf

.INTERMEDIATES: manuscripts/censored-mle.tex

models: stan/censored-mle-m0.rds \
	stan/censored-mle-m0-robust.rds \
	stan/censored-mle-m1.rds \
	stan/censored-mle-m2.rds \
	stan/censored-mle-m3.rds \
	stan/censored-mle-m4.rds

input-data: data/biofouling.rds data/imputations.rds

output-data: data/censored-mle-m0.rds \
	data/censored-mle-m0-robust.rds \
	data/censored-mle-m1.rds \
	data/censored-mle-m2.rds \
	data/censored-mle-m3.rds \
	data/censored-mle-m4.rds

PROC-DATA = graphics/obs-hist.pdf \
	graphics/imp-days1.pdf \
	graphics/imp-trips.pdf \
	graphics/imp-paint.pdf \
	graphics/plM1boat.pdf \
	graphics/plM1paint.pdf \
	graphics/plM3boat.pdf \
	graphics/plM3paint.pdf \
	graphics/plM4Type.pdf \
	graphics/plSummary.pdf \
	data/looic.rds

processed-data: $(PROC-DATA)

# Defaults for number of multiply imputed datasets and HMC iterations if not
# passed via cmdline.
NUMMI?=50
MCITER?=2000

################################################################################
# Make data for feeding into models and manuscript
data/imputations.rds data/biofouling.rds: scripts/data-cleaning.R \
	data-raw/samples.csv data-raw/vessel.csv
	cd $(<D); \
	Rscript --no-save --no-restore $(<F) numMI=$(NUMMI)

################################################################################
# Rules for making stan models
stan/censored-mle-m0.rds: scripts/compile-model.R stan/censored-mle-m0.stan
	cd $(<D); \
	Rscript --no-save --no-restore $(<F) mname=$(basename $(@F) .rds)

stan/censored-mle-m0-robust.rds: scripts/compile-model.R \
	stan/censored-mle-m0-robust.stan
	cd $(<D); \
	Rscript --no-save --no-restore $(<F) mname=$(basename $(@F) .rds)

stan/censored-mle-m1.rds: scripts/compile-model.R stan/censored-mle-m1.stan
	cd $(<D); \
	Rscript --no-save --no-restore $(<F) mname=$(basename $(@F) .rds)

stan/censored-mle-m2.rds: scripts/compile-model.R stan/censored-mle-m2.stan
	cd $(<D); \
	Rscript --no-save --no-restore $(<F) mname=$(basename $(@F) .rds)

stan/censored-mle-m3.rds: scripts/compile-model.R stan/censored-mle-m3.stan
	cd $(<D); \
	Rscript --no-save --no-restore $(<F) mname=$(basename $(@F) .rds)

stan/censored-mle-m4.rds: scripts/compile-model.R stan/censored-mle-m4.stan
	cd $(<D); \
	Rscript --no-save --no-restore $(<F) mname=$(basename $(@F) .rds)

################################################################################
# Rules to fit models with data
data/censored-mle-m0.rds: scripts/fit-model.R \
	stan/censored-mle-m0.rds data/imputations.rds
	cd $(<D); \
	Rscript --no-save --no-restore $(<F) mname=$(basename $(@F) .rds) \
		myseed=737 iter=$(MCITER)

data/censored-mle-m0-robust.rds: scripts/fit-model.R \
	stan/censored-mle-m0-robust.rds data/imputations.rds
	cd $(<D); \
	Rscript --no-save --no-restore $(<F) mname=$(basename $(@F) .rds) \
		myseed=737 iter=$(MCITER)

data/censored-mle-m1.rds: scripts/fit-model.R \
	stan/censored-mle-m1.rds data/imputations.rds
	cd $(<D); \
	Rscript --no-save --no-restore $(<F) mname=$(basename $(@F) .rds) \
		myseed=666 iter=$(MCITER)

data/censored-mle-m2.rds: scripts/fit-model.R \
	stan/censored-mle-m2.rds data/imputations.rds
	cd $(<D); \
	Rscript --no-save --no-restore $(<F) mname=$(basename $(@F) .rds) \
		myseed=42 iter=$(MCITER)

data/censored-mle-m3.rds: scripts/fit-model.R \
	stan/censored-mle-m3.rds data/imputations.rds
	cd $(<D); \
	Rscript --no-save --no-restore $(<F) mname=$(basename $(@F) .rds) \
		myseed=13 iter=$(MCITER)

data/censored-mle-m4.rds: scripts/fit-model.R \
	stan/censored-mle-m4.rds data/imputations.rds
	cd $(<D); \
	Rscript --no-save --no-restore $(<F) mname=$(basename $(@F) .rds) \
		myseed=987 iter=$(MCITER)

################################################################################
# Rules to process data (add dependencies later).
$(PROC-DATA): scripts/post-process.R output-data
	cd $(<D); \
	Rscript --no-save --no-restore $(<F)

################################################################################
# Rules to make manuscripts
manuscripts/censored-mle.html: manuscripts/censored-mle.Rmd \
	data-raw/samples.csv
	cd $(<D); \
	Rscript --no-save --no-restore -e "rmarkdown::render('$(<F)')"

manuscripts/model-interrogation.html: manuscripts/model-interrogation.Rmd \
	output-data
	cd $(<D); \
	Rscript --no-save --no-restore -e "rmarkdown::render('$(<F)')"

%.tex: %.Rnw data/biofouling.rds data/imputations.rds $(PROC-DATA)
	cd $(<D); \
	Rscript --no-save --no-restore -e "knitr::knit('$(<F)')"

%.pdf: %.tex
	cd $(<D); \
	latexmk -pdf $(<F)

################################################################################
# Cleaning targets
clean-models:
	cd stan/; \
	rm -f *.rds

clean-manuscripts:
	cd manuscripts/; \
	rm -rf *.aux *.bbl *.bcf *.blg *.fdb_latexmk *.fls *.lof *.log *.lot \
		*.code *.loe *.toc *.rec *.out *.run.xml *~ *.prv \
		censored-mle.tex _region_*

clobber: clean-data clean-manuscripts
	cd manuscripts/; \
	rm -rf auto/ cache/ figure/
