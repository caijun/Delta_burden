# system("g++ -v")
setwd("C:/Users/Jun Cai/Documents/C++/VOC_burden")

system("g++ svir_leaky_multiple.cpp -lgsl -o svir_leaky_multiple")


library(glue)
dout <- glue("./output/multiple/")

if(!dir.exists(dout)) {
  dir.create(dout, recursive = T)
}

# list.files(dout, full.names = T) -> tmp_files
# if (length(tmp_files) != 0) {
#   invisible(file.remove(tmp_files))
# }

# the vaccine efficacy after 14 days of the 2nd dose for wild-type and 4 variant strains; default is 80%
# ve2 <- c(0.8, 0.8, 0.8, 0.8, 0.8)
ve2 <- c(0.735, 0.667, 0.547, 0.493, 0.704)
write.table(ve2, file = "data/multiple/param/ve2", row.names = F, col.names = F)
# R0 for wild-type and 4 variant strains
R0 <- c(2.5, 4, 3.8, 5, 6)
write.table(R0, file = "data/multiple/param/R0", row.names = F, col.names = F)
# initial fractions of wild-type and 4 variant strains
# init_f <- c(0.2, 0.2, 0.2, 0.2, 0.2)
# only wild-type
# init_f <- c(1, 0, 0, 0, 0)
# only Delta
# init_f <- c(0, 0, 0, 0, 1)
# proportions of variants across globe
init_f <- c(0.113, 0.156, 0.005, 0.018, 0.708)
write.table(init_f, file = "data/multiple/param/init_f", row.names = F, col.names = F)

argv <- c()
argv[1:5] <- NA
argv[1] = dout # dirName
argv[2] = 275 # epidemic start day
argv[3] = 2 # vaccination strategy
argv[4] = 1 # capacity of daily vaccine doses
argv[5] = 200 # number of simulations (max nsim = 200)

t1 <- Sys.time()

exe.cmd <- paste("./svir_leaky_multiple", paste(argv, collapse = " "))
exe.cmd
system(exe.cmd)

t2 <- Sys.time()
t2 - t1
