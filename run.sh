#!/usr/bin/bash

R --no-save -e "rmarkdown::render('analysis.Rmd', 'html_document')"