library('googlesheets')
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library(doParallel))
cl = makeCluster(6)
registerDoParallel(cl, cores=6)

scratch='/scratch/groups/mschatz1/cpowgs/'
srcdir='~/Code/utils/marcc/'


gs_auth(token = "~/.ssh/googlesheets_token.rds")
status=gs_url('https://docs.google.com/spreadsheets/d/1K6rpcxvNPn0vX7qFKti150X9cYZtwPR-m_8TOC_RlTQ/edit#gid=0')
data=gs_read(status)

source('~/Code/kusama/update.R')
newdata=update_gsheet(data)

gs_edit_cells(status, input=newdata, anchor=paste0('A1'))




#Do all rows simultaneously: 
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
    if (!file.exists(paste0(datadir, '/untar_done.txt'))) {
        ##use appropriate untar script depending on how many there are
        if (is.na(data$raw2[i])) {
            system(paste0('sbatch --output=', batchlogsdir, '/untar_log.txt ~/Code/utils/marcc/untar.scr ', data$raw1[i], ' ', datadir))
        } else {
            system(paste0('sbatch --output=', batchlogsdir, '/untar_log.txt ~/Code/utils/marcc/untar2.scr ',  data$raw1[i], ' ', data$raw2[i], ' ', datadir))
        }
    }


    
    ##submit call jobs if untar is done, and calling hasn't been attempted yet
    numdirs=length(list.dirs(rawdir, recursive=FALSE))-1
    numcalled=length(list.files(calldonedir, recursive=FALSE))-1

    if (file.exists(paste0(datadir, '/untar_done.txt')) && numcalled!=numdirs) {
        system(paste0('mkdir -p ', calldir))
        system(paste0('mkdir -p ', calllogsdir))
        system(paste0('mkdir -p ', calldonedir))
        
        ##check how large the calling array needs to be, and submit
        system(paste0('sbatch --array=0-', as.character(numdirs), ' --job-name=', name, ' --output=', calllogsdir, '/', name, '.%A_%a.out ', srcdir , 'call.scr ', datadir))
    }



    ##submit fq jobs
    if (numcalled==numdirs && !file.exists(paste0(fqdir, '/', name, '.fq'))) {
        system(paste0('mkdir -p ', fqdir))
        
        ##submit the fq gather script
        system(paste0('sbatch ', '--output=',batchlogsdir, '/fq.out ', srcdir , 'fqs.scr ', datadir))

    }


    ##submit assembly job
    fq=paste0(fqdir, '/', name, '.fq')
    if (file.exists(fq)) {
        system(paste0('mkdir -p ', canudir))
        
        ##submit the fq gather script
        system(paste0('bash ',srcdir, '/assemble_bacteria.sh ', fq, ' ', canudir))
        
    }



}
