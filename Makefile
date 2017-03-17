# Time-stamp: <2017-03-17 10:05:45 (slane)>
.PHONY: all models input-data output-data clean-data clean-manuscripts clobber

all: manuscripts/censored-mle.html

# models: stan/dynamic-governance-m0.rds

# input-data: data/data-KIWI-10.rda

# output-data: data/stanfit-KIWI-dynamic-governance-m0-10.rda

################################################################################
# Rules for making stan models
# stan/dynamic-governance-m0.rds: R/compile-model.R
# 	cd $(<D); \
# 	Rscript $(<F) mname=dynamic-governance-m0

################################################################################
# Rules to make data for feeding into models
# data/data-KIWI-10.rda: R/create-data.R data/fruit.RData
# 	cd $(<D); \
# 	Rscript $(<F) pathway=KIWI size=10

################################################################################
# Rules to fit models with data
# data/stanfit-KIWI-dynamic-governance-m0-10.rda: R/fit-model.R \
# 	stan/dynamic-governance-m0.rds data/data-KIWI-10.rda
# 	cd $(<D); \
# 	Rscript $(<F) mname=dynamic-governance-m0 pathway=KIWI size=10 iter=2000

################################################################################
# Rules to make manuscripts
%.html: %.Rmd data/biofouling.csv
	cd $(<D); \
	echo "rmarkdown::render('$(<F)')" | R --no-save --no-restore

# %.tex: %.Rnw data/stanfit-KIWI-dynamic-governance-m1-10.rda
# 	cd $(<D); \
# 	echo "knitr::knit('$(<F)')" | R --no-save --no-restore

# %.pdf: %.tex
# 	cd $(<D); \
# 	latexmk -pdf $(<F)

################################################################################
# Cleaning targets
clean-data:
	cd stan/; \
	rm -f *.rds

clean-manuscripts:
	cd manuscripts/; \
	rm -f *.aux *.bbl *.bcf *.blg *.fdb_latexmk *.fls *.lof *.log *.lot \
		*.code *.loe *.toc *.rec *.out *.run.xml *~

clobber:
	cd manuscripts/; \
	rm -rf *.aux *.bbl *.bcf *.blg *.fdb_latexmk *.fls *.lof *.log *.lot \
		*.code *.loe *.toc *.rec *.out *.run.xml *~ auto/ cache/ figure/
