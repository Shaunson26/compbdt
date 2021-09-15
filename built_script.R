library(devtools)

#use_readme_rmd()

load_all()
document()

# Plan
# Split up compdbt into components, which are?
# input checking
# function body
# output

# R/compbdt_full
# comment out parts that have been done

# Function files ----
#use_r('compbdt_full')

# Check any or all zero
#use_r('check_zero')
check_zero(c(0,1), 'any')
check_zero(c(0,0), 'all')
# Check values less than zero
