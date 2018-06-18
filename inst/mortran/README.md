## Generating the fortran file

In an R session, source the `process.R` file. 

This will 
- Replace the reals with double precision
- Notify about any long lines, just in case
- Run mortran processor
- Chop the lines at 72 chars 
- Finally produce glasso.f
