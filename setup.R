library(packrat)
# Set a working directory to the application path: ~/ShinyApps/<R_VERSION>/<APP_NAME>
setwd('~/ShinyApps/4.4.1/ngs_metadata/')

# Initialize 'packrat' directory structure and environment
packrat::init(options = list(symlink.system.packages = FALSE))

# Embed application dependencies. Update a list in 'app.deps' variable to address the specific application needs
app.deps <- c("shiny", "shinydashboard","readxl","rhandsontable","dplyr","packrat")

install.packages(app.deps, repos="http://cran.rstudio.com/", dependencies=TRUE)
