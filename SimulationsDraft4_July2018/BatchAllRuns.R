#devtools::install_github("bomeara/BMhyb")

n.cores=parallel::detectCores()
for (core.index in sequence(n.cores)) {
	system("nohup R CMD BATCH BatchRuns.R > /dev/null &")
	Sys.sleep(5)
}
