# system("g++ svir_AON_single.cpp -lgsl -o svir_AON_single")

#------ load packages ------
library(rstudioapi)
library(glue)
Sys.setlocale("LC_TIME", "US")

#------ generate scripts ------
outdir <- glue("output/single/AON/batch/main/")

if(!dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}
list.files(outdir, full.names = T) -> tmp_files
if (length(tmp_files) != 0) {
  invisible(file.remove(tmp_files))
}

nsim <- 200

# main analysis ################################################################
scripts <- c()
R0 <- 6
ve2 <- 0.543
for (Tonset in c(275, 305, 336)) {
  for (capacity in c(1)) {
    for (strategy in 1:2) {
      for (susflag in 1:2) {
        for (cmflag in 1:2) {
          if (!(susflag == 2 && cmflag == 2)) {
            exe.cmd <- paste("svir_AON_single.exe", outdir, Tonset, strategy, capacity, ve2, R0, nsim, susflag, cmflag)
            scripts <- c(scripts, exe.cmd)
          }
        }
      }
    }
  }
}


# # increase intensity of NPIs and VE ############################################
# scripts <- c()
# R0max <- 6.1
# for (R0 in seq(1.1, R0max, by = 0.2)) {
#   for (ve2 in seq(0.5, 0.95, by = 0.05)) {
#     for (Tonset in c(275, 305)) {
#       for (capacity in c(1)) {
#         for (strategy in 1:2) {
#           for (susflag in c(1)) {
#             for (cmflag in c(1)) {
#               exe.cmd <- paste("svir_AON_single.exe", outdir, Tonset, strategy, capacity, ve2, R0, nsim, susflag, cmflag)
#               scripts <- c(scripts, exe.cmd)
#             }
#           }
#         }
#       }
#     }
#   }
# }


#------ batch run ------
t1 <- Sys.time()

max_run <- 8
sleep_time <- 2

if (length(scripts) < max_run) {
  max_run <- length(scripts)
}

log_save <- c()
terminal_id <- c()
for (i in 1:length(scripts)) {
  terminalExecute(scripts[i], show = F) -> term_id
  c(terminal_id, term_id) -> terminal_id
  data.frame(i = i, script = scripts[i], terminal_id = term_id, finish = F) -> log_tmp
  rbind(log_save, log_tmp) -> log_save
  
  length(terminal_id) -> tmp_run
  if (tmp_run < max_run & i <= length(scripts)) {
    next
  }
  
  num_exit <- 0
  while (num_exit == 0) {
    terminal_exited <- c()
    for (j in terminal_id) {
      c(terminal_exited, !is.null(terminalExitCode(j))) -> terminal_exited
    }
    num_exit <- sum(terminal_exited)
    Sys.sleep(sleep_time)
  }
  
  terminal_id[terminal_exited] -> terminal_need_2_kill
  
  for (k in terminal_need_2_kill) {
    terminalKill(k)
  }
  terminal_id[!terminal_exited] -> terminal_id
  
  if(i == length(scripts)) {
    length(terminal_id) -> tmp_run
    num_exit <- 0
    while (num_exit != tmp_run) {
      terminal_exited <- c()
      for (j in terminal_id) {
        c(terminal_exited, !is.null(terminalExitCode(j))) -> terminal_exited
      }
      num_exit <- sum(terminal_exited)
      Sys.sleep(sleep_time)
    }
    
    terminal_id[terminal_exited] -> terminal_need_2_kill
    
    for (k in terminal_need_2_kill) {
      terminalKill(k)
    }
    terminal_id[!terminal_exited] -> terminal_id
  }
  log_save$finish[log_save$terminal_id %in% terminal_need_2_kill] <- T
}

t2 <- Sys.time()
t2 - t1
