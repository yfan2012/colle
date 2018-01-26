library('googlesheets')
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library(doParallel))

gs_auth(token = "~/.ssh/googlesheets_token.rds")
status=gs_url('https://docs.google.com/spreadsheets/d/1BbVvmVUaJGPFUkUO22Go9VMrLZ9iOojeND__E2SD4LY/edit#gid=0')

data=gs_read(status)

cl = makeCluster(6)
registerDoParallel(cl, cores=6)

scratch='/scratch/groups/mschatz1/cpowgs/'
srcdir='~/Code/carbapenem_r21/auto_marcc/'



##Do all rows simultaneously: 
foreach(i=1:dim(data)[1]) %dopar% {
    name=paste0(data$org[i], '_' , data$isolate_num[i])

    ##untar raw
    if (!dir.exists(paste0(scratch, name))) {
        system(paste0('mkdir -p ', scratch, name))
        system(paste0('mkdir -p ', scratch, name, '/raw'))
        ##combine minion tars if necessary
        if (data$raw2[i]!=NA) {
            system(paste0('tar --concatenate file=', raw1, ' ', raw2))
        }
        system(paste0('sbatch ', srcdir, untar.scr, ' ', scratch, name, ' ', raw1))
    }
    
}
