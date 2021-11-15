# system("g++ -v")
setwd("C:/Users/Jun Cai/Documents/C++/VOC_burden")

system("g++ svir_leaky_single.cpp -lgsl -o svir_leaky_single")


library(glue)
dout <- glue("./output/single/leaky/")

if(!dir.exists(dout)) {
  dir.create(dout, recursive = T)
}

# list.files(dout, full.names = T) -> tmp_files
# if (length(tmp_files) != 0) {
#   invisible(file.remove(tmp_files))
# }

argv <- c()
argv[1:11] <- NA
argv[1] = dout # dirName
argv[2] = 366 # epidemic start day (2021-12-01)
argv[3] = 2 # vaccination strategy
argv[4] = 1 # capacity of daily vaccine doses
argv[5] = 0.543 # vaccine efficacy after 14 days of the 2nd dose
argv[6] = 6 # R0
argv[7] = 200 # number of simulations (max nsim = 200)
argv[8] = 1 # susceptibility to infection by age (1 = heterogeneous, 2 = homogeneous)
argv[9] = 1 # in which period contact matrix will be used (1 = baseline, 2 = postlockdown)
argv[10] = 40 # number of initial seed infectors
argv[11] = 7 # generation time (days)

t1 <- Sys.time()

exe.cmd <- paste("./svir_leaky_single", paste(argv, collapse = " "))
exe.cmd
system(exe.cmd)

t2 <- Sys.time()
t2 - t1
