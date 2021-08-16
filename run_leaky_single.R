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
argv[1:9] <- NA
argv[1] = dout # dirName
argv[2] = 275 # epidemic start day
argv[3] = 2 # vaccination strategy
argv[4] = 1 # capacity of daily vaccine doses
argv[5] = 0.8 # vaccine efficacy after 14 days of the 2nd dose
argv[6] = 6 # R0
argv[7] = 200 # number of simulations (max nsim = 200)
argv[8] = 2 # susceptibility to infection by age (1 = heterogeneous, 2 = homogeneous)
argv[9] = 2 # in which period contact matrix will be used (1 = baseline, 2 = postlockdown)

t1 <- Sys.time()

exe.cmd <- paste("./svir_leaky_single", paste(argv, collapse = " "))
exe.cmd
system(exe.cmd)

t2 <- Sys.time()
t2 - t1
