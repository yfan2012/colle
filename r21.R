library('googlesheets')
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(plyr))
cl = makeCluster(6)
registerDoParallel(cl, cores=6)

work='/work-zfs/mschatz1/cpowgs/r21/'
scratch='/scratch/groups/mschatz1/cpowgs/r21/'
srcdir='~/Code/utils/marcc/'

gs_auth(token = "~/.ssh/googlesheets_token.rds")
status=gs_url('https://docs.google.com/spreadsheets/d/1BbVvmVUaJGPFUkUO22Go9VMrLZ9iOojeND__E2SD4LY/edit#gid=0')
data=gs_read(status)

source('~/Code/colle/update.R')
newdata=update_gsheet(data, scratch)

gs_edit_cells(status, input=newdata, anchor=paste0('A1'))

data=newdata


##Do all rows simultaneously: 
##foreach(i=1:dim(data)[1]) %dopar% {


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
    bamdir=paste0(datadir, '/bams')
    canudir=paste0(datadir, '/canu_assembly')
    polishdir=paste0(datadir, '/polish')
    mpolishdir=paste0(datadir, '/mpolish')
    pilondir=paste0(datadir,'/pilon')
    
    system(paste0('mkdir -p ', datadir))
    system(paste0('mkdir -p ', batchlogsdir))
    system(paste0('mkdir -p ', rawdir))
    
    
    ##untar raw if raw dir isn't there
    if (!file.exists(paste0(datadir, '/untar_done.txt')) && !identical(x$untar, 'submitted')) {
        ##use appropriate untar script depending on how many there are
        if (is.na(x$raw2)) {
            system(paste0('sbatch --output=', batchlogsdir, '/untar_log.txt ~/Code/utils/marcc/untar.scr ', paste0(work,x$raw1), ' ', datadir))
            x$untar='submitted'
        } else {
            system(paste0('sbatch --output=', batchlogsdir, '/untar_log.txt ~/Code/utils/marcc/untar2.scr ',  paste0(work, x$raw1), ' ', paste0(work,x$raw2), ' ', datadir))
            x$untar='submitted'
        }
    }
    
    
    
    ##submit call jobs if untar is done, and calling hasn't been attempted yet
    numdirs=length(list.dirs(rawdir, recursive=FALSE))-1
    numcalled=length(list.files(calldonedir, recursive=FALSE))-1
    fq=paste0(fqdir, '/', name, '.fq')    
    if (file.exists(paste0(datadir, '/untar_done.txt')) && numcalled!=numdirs && !identical(x$call, 'submitted') && !file.exists(fq) ) {
        system(paste0('mkdir -p ', calldir))
        system(paste0('mkdir -p ', calllogsdir))
        system(paste0('mkdir -p ', calldonedir))
        
        ##check how large the calling array needs to be, and submit
        system(paste0('sbatch --array=0-', as.character(numdirs), ' --job-name=', name, ' --output=', calllogsdir, '/', name, '.%A_%a.out ', srcdir , 'call.scr ', datadir))
        x$call='submitted'
    }
    
    
    ##submit fq jobs
    if (numcalled==numdirs && numdirs!=-1 && !file.exists(fq) && !identical(x$fq, 'submitted')) {
        system(paste0('mkdir -p ', fqdir))
        
        ##submit the fq gather script
        system(paste0('sbatch ', '--output=',batchlogsdir, '/fq.out ', srcdir , 'fqs.scr ', datadir))
        x$fq='submitted'
        
    }

    
    ##submit assembly job
    assembly=paste0(canudir, '/', name, '.contigs.fasta')
    if (file.exists(fq) && !identical(x$assembly, 'submitted') && !file.exists(assembly)) {
        system(paste0('mkdir -p ', canudir))
        
        ##submit the assembly script
        system(paste0('bash ',srcdir, '/assemble_bacteria.sh ', fq, ' ', canudir))
        x$assembly='submitted'
    }

    
    ##alignment
    alignment=paste0(bamdir, '/', name, '.sorted.bam')
    if (file.exists(fq) && !identical(x$align, 'submitted') && file.exists(assembly) && !file.exists(alignment)) {
        system(paste0('mkdir -p ', bamdir))

        ##submit the align script
        system(paste0('sbatch ', '--output=', batchlogsdir, '/align.out ', srcdir, 'align.scr ', datadir, ' ', assembly))
        x$align='submitted'
    }
    

    ##polish
    polished=paste0(polishdir, '/', name, '.polished.fasta')
    if (file.exists(fq) && file.exists(assembly) && file.exists(alignment) && !identical(x$polish, 'submitted')) {
        system(paste0('mkdir -p ', polishdir))
        
        if (!file.exists(polished) || file.info(polished)$size==0) {
            ##submit the polish script
            system(paste0('sbatch --output=', batchlogsdir, '/polish.out --job-name=', name,' ', srcdir, 'polish.scr ', datadir))
            x$polish='submitted'
        }
    }

    
    ##mpolish
    mpolished=paste0(mpolishdir, '/', name, '.polished_meth.fasta')
    if (file.exists(fq) && file.exists(assembly) && file.exists(alignment) && !identical(x$mpolish, 'submitted') && file.exists(polished)) {
        system(paste0('mkdir -p ', mpolishdir))

        if (!file.exists(mpolished) || (file.exists(mpolished) && file.info(mpolished)$size==0)) {
            ##submit the polish script
            system(paste0('sbatch --output=', batchlogsdir, '/polish_meth.out --job-name=', name,' ', srcdir, 'polish_meth.scr ', datadir))
            x$mpolish='submitted'
        }
    }
    

    ##pilon
    pilon=paste0(pilondir, '/', name, '.pilon.10.fasta')
    illpre=paste0(x$org, '-' , x$isolate_num)

    
    if (file.exists(assembly) && !file.exists(pilon) && !identical(x$pilon, 'submitted')) {
        ##submit the pilon script - script assumes location of illumina data
        system(paste0('rm -rf ', pilondir, '/*'))
        system(paste0('mkdir -p ', pilondir))
        system(paste0('cp /work-zfs/mschatz1/cpowgs/illumina/klpn_all/', illpre, '_*.fastq.gz ', pilondir, '/'))
        system(paste0('rename - _ ', pilondir, '/*'))
        illreads=list.files(pilondir, pattern='fastq.gz')

        if (length(illreads)>1){
            system(paste0('sbatch --output=', batchlogsdir, '/pilon.out --job-name=', name, ' ', srcdir, 'pilon.scr ', datadir)) 
            x$pilon='submitted'
        }
    }
    
    
    return(x)        
})

gs_edit_cells(status, input=updata, anchor=paste0('A1'))
