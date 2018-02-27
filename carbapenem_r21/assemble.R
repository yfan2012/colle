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


for (i in 1:dim(data)[1]){
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

    ##check if untar is done
    if (file.exists(paste0(datadir, '/untar_done.txt'))) {
        data$untar[i]='done'
    }

    ##check if calling isdone
    numdirs=length(list.dirs(rawdir, recursive=FALSE))-1
    numcalled=length(list.files(calldonedir, recursive=FALSE))-1
    if (numdirs==numcalled) {
        data$call[i]='done'
    }

    ##check if fq exits
    if (file.exists(paste0(fqdir,'/', name, '.fq'))) {
        data$fq[i]='done'
    }

    
    if (file.exists(paste0(canudir,'/', name, '.contigs.fasta'))) {
        data$assembly[i]='done'
    }
}
gs_edit_cells(status, input=data, anchor=paste0('A1'))




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

    system(paste0('mkdir -p ', canudir))

    numdirs=length(list.dirs(rawdir, recursive=FALSE))-1
    numcalled=length(list.files(calldonedir, recursive=FALSE))-1

    ##submit assembly of only passing reads
    fq=paste0(fqdir, '/', name, '.fq')
    if (file.exists(fq)) {

        ##submit the fq gather script
        system(paste0('bash ',srcdir, '/assemble_bacteria.sh ', fq, ' ', canudir))

    }
}
