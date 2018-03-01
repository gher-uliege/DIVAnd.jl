figdir = "./figures/"
outputdir = "./netCDF/"
logging
isdir(figdir) ? info("Directory already exists") : mkdir(figdir);
isdir(outputdir) ? info("Directory already exists") : mkdir(outputdir);
