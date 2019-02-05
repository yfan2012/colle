library('googlesheets')
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(plyr))
cl = makeCluster(6)
registerDoParallel(cl, cores=6)

work='/work-zfs/mschatz1/cpowgs/'
scratch='/scratch/groups/mschatz1/cpowgs/'
srcdir='~/Code/utils/marcc/'

gs_auth(token = "~/.ssh/googlesheets_token.rds")
status=gs_url('https://docs.google.com/spreadsheets/d/1BbVvmVUaJGPFUkUO22Go9VMrLZ9iOojeND__E2SD4LY/edit#gid=0')
data=gs_read(status)

source('~/Code/colle/update.R')
newdata=update_gsheet(data, scratch)

gs_edit_cells(status, input=newdata, anchor=paste0('A1'))

data=newdata

##foreach(i=sub) %dopar% {
updata=ddply(data, .(isolate_num), function(x) {
    name=paste0(x$org, '_' , x$isolate_num)
    
    datadir=paste0(scratch, name)
    rawdir=paste0(datadir,'/raw')
    calldir=paste0(datadir,'/called')
    calllogsdir=paste0(datadir,'/call_logs')
    calldonedir=paste0(datadir,'/call_done')
    batchlogsdir=paste0(datadir,'/batch_logs')
    fqdir=paste0(datadir,'/fastqs')
    bamdir=paste0(datadir, '/bams17')
    canudir=paste0(datadir, '/canu17')
    polishdir=paste0(datadir, '/polish17')
    mpolishdir=paste0(datadir, '/mpolish17')
    pilondir=paste0(datadir,'/pilon17')
    
    if (identical(x$call,'done') && identical(x$polish,'done') && identical(x$mpolish,'done')) {
        system(paste0('rm -r ', rawdir))
    }
})

