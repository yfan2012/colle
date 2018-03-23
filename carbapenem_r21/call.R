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

source('~/Code/kusama/update.R')
newdata=update_gsheet(data)

gs_edit_cells(status, input=newdata, anchor=paste0('A1'))




##Call all rows simultaneously: 
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

    system(paste0('mkdir -p ', calldir))
    system(paste0('mkdir -p ', calllogsdir))
    system(paste0('mkdir -p ', calldonedir))

    numdirs=length(list.dirs(rawdir, recursive=FALSE))-1
    numcalled=length(list.files(calldonedir, recursive=FALSE))-1


    ##submit call jobs if untar is done, and calling hasn't been attempted yet
    if (file.exists(paste0(datadir, '/untar_done.txt')) && numcalled!=numdirs) {

        system(paste0('mkdir -p ', calldir))
        system(paste0('mkdir -p ', calllogsdir))
        system(paste0('mkdir -p ', calldonedir))
        
        ##check how large the calling array needs to be, and submit
        system(paste0('sbatch --array=0-', as.character(numdirs), ' --job-name=', name, ' --output=', calllogsdir, '/', name, '.%A_%a.out ', srcdir , 'call.scr ', datadir))

    }
}
