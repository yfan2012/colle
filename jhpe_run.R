library('googlesheets')
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library(doParallel))
cl = makeCluster(6)
registerDoParallel(cl, cores=6)

epidir='/work-zfs/mschatz1/cpowgs/epicenter/'
scratch='/scratch/groups/mschatz1/cpowgs/epicenter/'
srcdir='~/Code/utils/marcc/'

gs_auth(token = "~/.ssh/googlesheets_token.rds")
status=gs_url('https://docs.google.com/spreadsheets/d/1K6rpcxvNPn0vX7qFKti150X9cYZtwPR-m_8TOC_RlTQ/edit#gid=0')
olddata=gs_read(status)
olddata$barcode=gsub('NB', 'barcode', olddata$barcode)

##check to see if stuff has finished running since last time
source('~/Code/colle/jhpe_update.R')
data=update_gsheet(olddata, scratch)
gs_edit_cells(status, input=data, anchor=paste0('A1'))

allruns=unique(data$run_code)

##foreach(i=allruns) %dopar% {
calldata=ddply(data, .(run_code), function(x) {
    
    rundir=paste0(scratch, x$run_name[1])
    rundir_out=paste0(rundir, '/batch_logs')
    rundir_calllogs=paste0(rundir, '/call_logs')
    rundir_raw=paste0(rundir, '/raw')
    rundir_call=paste0(rundir, '/called')
    rundir_done=paste0(rundir, '/call_done')

    system(paste0('mkdir -p ', rundir))
    system(paste0('mkdir -p ', rundir_out))
    system(paste0('mkdir -p ', rundir_raw))
    
    ##untar
    untar_done=paste0(rundir, '/untar_done.txt')
    if (!file.exists(untar_done) && !identical(x$untar[1],'submitted') && !is.na(x$raw1[1])) {
        if (is.na(x$raw2[1])) {
            system(paste0('sbatch --output=', rundir_out, '/untar_log.txt ', srcdir, '/untar.scr ', epidir, x$raw1[1], ' ', rundir))
        } else {
            system(paste0('sbatch --output=', rundir_out, '/untar_log.txt ', srcdir, '/untar2.scr ', epidir, x$raw1[1], ' ', epidir, x$raw2[1], ' ', rundir))
        }
        x$untar='submitted'
    }
    

    
    ##submit calling jobs if untar is done
    numdirs=length(list.dirs(rundir_raw, recursive=FALSE))-1
    numcalled=length(list.files(rundir_done, recursive=FALSE))-1
    
    if (file.exists(untar_done) && numcalled!=numdirs && !identical(x$call[1],'submitted')) { 

        system(paste0('mkdir -p ', rundir_call))
        system(paste0('mkdir -p ', rundir_calllogs))
        system(paste0('mkdir -p ', rundir_done))

        system(paste0('sbatch --array=0-', as.character(numdirs), ' --job-name=', x$run_name[1], ' --output=', rundir_calllogs, '/', x$run_name[1], '.%A_%a.out ', srcdir, 'bc_call.scr ', rundir))
        
        x$call='submitted'
    }

    return(x)
})


updata=ddply(calldata, 1, function(x){
    name=paste0(x$org,'_', x$isolate_num)
    
    datadir=paste0(scratch,name)
    rundir=paste0(scratch, x$run_name)
    
    rawdir=paste0(rundir, '/raw')
    calldir=paste0(rundir, '/called')
    calldonedir=paste0(rundir, '/call_done')
    batchlogsdir=paste0(datadir,'/batch_logs')
    fqdir=paste0(datadir, '/fastqs')
    canudir=paste0(datadir,'/canu_assembly')    
    bamdir=paste0(datadir, '/bams')
    polishdir=paste0(datadir, '/polish')
    mpolishdir=paste0(datadir, '/mpolish')
    

    system(paste0('mkdir -p ', datadir))
    system(paste0('mkdir -p ', batchlogsdir))

    numdirs=length(list.dirs(rawdir, recursive=FALSE))-1
    numcalled=length(list.files(calldonedir, recursive=FALSE))-1
    
    ##gather all the fqs if calling is done
    fq=paste0(fqdir, '/', name, '.fq')        
    if (numdirs==numcalled && numdirs!=-1 && !file.exists(fq) && !identical(x$fq, 'submitted')) {
        system(paste0('mkdir -p ', fqdir))
        barcode=x$barcode
        system(paste0('sbatch --output=', batchlogsdir, '/fq.out ', srcdir, 'bc_fqs.scr ', rundir, ' ', fqdir, ' ', barcode))
        x$fq='submitted'
    }


    ##submit assembly
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
        system(paste0('sbatch ', '--output=', batchlogsdir, '/align.out ', srcdir, 'align_old.scr ', datadir, ' ', assembly))
        x$align='submitted'
    }
    
    
    ##polish
    polished=paste0(polishdir, '/', name, '.polished.fasta')
    if (file.exists(fq) && file.exists(assembly) && file.exists(alignment) && !identical(x$polish, 'submitted')) {
        system(paste0('mkdir -p ', polishdir))
        
        ##submit the polish script only if polished does not exist or if polished is size0
        if (!file.exists(polished) || file.info(polished)$size==0) {
            system(paste0('sbatch --output=', batchlogsdir, '/polish.out --job-name=', name,' ', srcdir, 'polish.scr ', datadir, ' ', rundir))
            x$polish='submitted'
        }
    }
    
    
    ##mpolish
    mpolished=paste0(mpolishdir, '/', name, '.polished_meth.fasta')
    if (file.exists(fq) && file.exists(assembly) && file.exists(alignment) && !identical(x$mpolish, 'submitted') && file.exists(polished)) {
        system(paste0('mkdir -p ', mpolishdir))
        
        ##submit the polish script
        if (!file.exists(mpolished) || file.info(mpolished)$size==0) {
            system(paste0('sbatch --output=', batchlogsdir, '/polish_meth.out --job-name=', name,' ', srcdir, 'polish_meth.scr ', datadir))
            x$mpolish='submitted'
        }
    }
    
    
    
    
    return(x)
})



gs_edit_cells(status, input=updata, anchor=paste0('A1'))
