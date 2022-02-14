# system("g++ svir_leaky_single.cpp -lgsl -o svir_leaky_single")

#------ load packages ------
library(rstudioapi)
library(glue)
library(tidyverse)
Sys.setlocale("LC_TIME", "US")

#------ generate scripts ------
outdir <- glue("output/single/leaky/batch/sup3/")

if(!dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}
list.files(outdir, full.names = T) -> tmp_files
if (length(tmp_files) != 0) {
  invisible(file.remove(tmp_files))
}

nsim <- 200

# # main analysis ################################################################
# scripts <- c()
# experiments <- expand.grid(R0 = 6, ve2 = 0.543, Tonset = 366, strategy = 1:2,  capacity = 1,
#                            susflag = 1:2, cmflag = 1:2, ni = c(40, 20, 10), gt = c(7.0, 4.6)) %>%
#   filter((cmflag == 1 & ni == 40 & gt == 7.0) | (susflag == 1 & ni == 40 & gt == 7.0) |
#            (susflag == 1 & cmflag == 1 & gt == 7.0) | (susflag == 1 & cmflag == 1 & ni == 40))
# 
# for(i in 1:nrow(experiments)) {
#   Tonset <- experiments$Tonset[i]
#   strategy <- experiments$strategy[i]
#   capacity <- experiments$capacity[i]
#   ve2 <- experiments$ve2[i]
#   R0 <- experiments$R0[i]
#   susflag <- experiments$susflag[i]
#   cmflag <- experiments$cmflag[i]
#   ni <- experiments$ni[i]
#   gt <- experiments$gt[i]
#   exe.cmd <- paste("svir_leaky_single.exe", outdir, Tonset, strategy, capacity, ve2, R0, nsim, susflag, cmflag, ni, gt)
#   scripts <- c(scripts, exe.cmd)
# }


# # no vax ################################################################
# scripts <- c()
# experiments <- expand.grid(R0 = 6, ve2 = 0, Tonset = 366, strategy = 1:2,  capacity = 1,
#                            susflag = 1:2, cmflag = 1:2, ni = c(40, 20, 10), gt = c(7.0, 4.6)) %>%
#   filter((cmflag == 1 & ni == 40 & gt == 7.0) | (susflag == 1 & ni == 40 & gt == 7.0) |
#            (susflag == 1 & cmflag == 1 & gt == 7.0) | (susflag == 1 & cmflag == 1 & ni == 40))
# 
# for(i in 1:nrow(experiments)) {
#   Tonset <- experiments$Tonset[i]
#   strategy <- experiments$strategy[i]
#   capacity <- experiments$capacity[i]
#   ve2 <- experiments$ve2[i]
#   R0 <- experiments$R0[i]
#   susflag <- experiments$susflag[i]
#   cmflag <- experiments$cmflag[i]
#   ni <- experiments$ni[i]
#   gt <- experiments$gt[i]
#   exe.cmd <- paste("svir_leaky_single.exe", outdir, Tonset, strategy, capacity, ve2, R0, nsim, susflag, cmflag, ni, gt)
#   scripts <- c(scripts, exe.cmd)
# }


# # increase intensity of NPIs and VE ############################################
# scripts <- c()
# R0max <- 6.1
# for (R0 in seq(1.1, R0max, by = 0.2)) {
#   for (ve2 in seq(0.5, 0.95, by = 0.05)) {
#     for (Tonset in c(366)) {
#       for (capacity in c(1)) {
#         for (strategy in c(2)) {
#           for (susflag in c(1)) {
#             for (cmflag in c(1)) {
#               for (ni in c(40)) {
#                 for (gt in c(7.0)) {
#                   exe.cmd <- paste("svir_leaky_single.exe", outdir, Tonset, strategy, capacity, ve2, R0, nsim, susflag, cmflag, ni, gt)
#                   scripts <- c(scripts, exe.cmd)
#                 }
#               }
#             }
#           }
#         }
#       }
#     }
#   }
# }


# # supplementary analysis1 ############################################
# scripts <- c()
# for (R0 in c(4.0, 4.1, 4.2, 4.3, 4.4, 4.5)) {
#   for (ve2 in c(0.75)) {
#     for (Tonset in c(366)) {
#       for (capacity in c(1)) {
#         for (strategy in c(2)) {
#           for (susflag in c(1)) {
#             for (cmflag in c(1)) {
#               for (ni in c(40)) {
#                 for (gt in c(7.0)) {
#                   exe.cmd <- paste("svir_leaky_single.exe", outdir, Tonset, strategy, capacity, ve2, R0, nsim, susflag, cmflag, ni, gt)
#                   scripts <- c(scripts, exe.cmd)
#                 }
#               }
#             }
#           }
#         }
#       }
#     }
#   }
# }


# # supplementary analysis2 ############################################
# scripts <- c()
# for (R0 in c(2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7)) {
#   for (ve2 in c(0.543)) {
#     for (Tonset in c(366)) {
#       for (capacity in c(1)) {
#         for (strategy in c(2)) {
#           for (susflag in c(1)) {
#             for (cmflag in c(1)) {
#               for (ni in c(40)) {
#                 for (gt in c(7.0)) {
#                   exe.cmd <- paste("svir_leaky_single.exe", outdir, Tonset, strategy, capacity, ve2, R0, nsim, susflag, cmflag, ni, gt)
#                   scripts <- c(scripts, exe.cmd)
#                 }
#               }
#             }
#           }
#         }
#       }
#     }
#   }
# }


# supplementary analysis3 ############################################
scripts <- c()
for (R0 in c(6.0)) {
  for (ve2 in c(0.80, 0.81, 0.82, 0.83, 0.84, 0.85)) {
    for (Tonset in c(366)) {
      for (capacity in c(1)) {
        for (strategy in c(2)) {
          for (susflag in c(1)) {
            for (cmflag in c(1)) {
              for (ni in c(40)) {
                for (gt in c(7.0)) {
                  exe.cmd <- paste("svir_leaky_single.exe", outdir, Tonset, strategy, capacity, ve2, R0, nsim, susflag, cmflag, ni, gt)
                  scripts <- c(scripts, exe.cmd)
                }
              }
            }
          }
        }
      }
    }
  }
}


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
