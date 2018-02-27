library('googlesheets')
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library(doParallel))


cl = makeCluster(6)
registerDoParallel(cl, cores=6)

scratch='/scratch/groups/mschatz1/cpowgs/'
srcdir='~/Code/utils/marcc/'

gs_auth(token = "~/.ssh/googlesheets_token.rds")
status=gs_url('https://docs.google.com/spreadsheets/d/1BbVvmVUaJGPFUkUO22Go9VMrLZ9iOojeND__E2SD4LY/edit#gid=0')
data=gs_read(status)


##Do all rows simultaneously: 
foreach(i=1:dim(data)[1]) %dopar% {
    name=paste0(data$org[i], '_' , data$isolate_num[i])

    datadir=paste0(scratch,name)
    rawdir=paste0(datadir,'/raw')
    calldir=paste0(datadir,'/called')
    calllogsdir=paste0(datadir,'/call_logs')
    calldonedir=paste0(datadir,'/call_done')
    batchlogsdir=paste0(datadir,'/batch_logs')
    fqdir=paste0(datadir,'/fastqs')
    canudir=paste0(datadir, '/canu_assembly')
    polishdir=paste0(datadir, '/polish')
    

    system(paste0('mkdir -p ', datadir))
    system(paste0('mkdir -p ', batchlogsdir))
    system(paste0('mkdir -p ', rawdir))
    
    ##untar raw if raw dir isn't there
    if (!file.exists(paste0(datadir, 'untar_done.txt'))) {

        
        ##use appropriate untar script depending on how many there are
        if (is.na(data$raw2[i])) {
            system(paste0('sbatch --output=', batchlogsdir, '/untar_log.txt ~/Code/utils/marcc/untar.scr ', data$raw1[i], ' ', datadir))
        } else {
            system(paste0('sbatch --output=', batchlogsdir, '/untar_log.txt ~/Code/utils/marcc/untar2.scr ',  data$raw1[i], ' ', data$raw2[i], ' ', datadir))
        }
    }

}
