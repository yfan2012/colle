library('googlesheets')
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library(doParallel))
cl = makeCluster(6)
registerDoParallel(cl, cores=6)

epidir='/work-zfs/mschatz1/cpowgs/epicenter/'
datadir='/scratch/groups/mschatz1/cpowgs/epicenter/'
srcdir='~/Code/utils/marcc/'

gs_auth(token = "~/.ssh/googlesheets_token.rds")
status=gs_url('https://docs.google.com/spreadsheets/d/1K6rpcxvNPn0vX7qFKti150X9cYZtwPR-m_8TOC_RlTQ/edit#gid=0')
olddata=gs_read(status)


##check to see if stuff has finished running since last time
source('~/Code/kusama/update.R')
data=update_gsheet(olddata, datadir)
gs_edit_cells(status, input=data, anchor=paste0('A1'))

allruns=unique(data$run_code)



foreach(i=allruns) %dopar% {
    runsamps=data[data$run_code==i,]
    
    rundir=paste0(datadir, runsamps$run_name[1])
    rundir_out=paste0(rundir, '/batch_logs')
    rundir_raw=paste0(rundir, '/raw')
    rundir_call=paste0(rundir, '/called')
    rundir_done=paste0(rundir, '/call_done' )
    
    system(paste0('mkdir -p ', rundir))
    system(paste0('mkdir -p ', rundir_out))
    system(paste0('mkdir -p ', rundir_raw))


    


    ##untar
    untar_done=paste0(rundir, '/untar_done.txt')
    if (!file.exists(untar_done) && runsamps$untar[1]!='submitted') {
        if (is.na(runsamps$raw2[1])) {
            system(paste0('sbatch --output=', rundir_out, '/untar_log.txt ', srcdir, '/untar.scr ', runsamps$raw1[1], ' ', rundir))
        } else {
            system(paste0('sbatch --output= ', rundir_out, '/untar_log.txt ', srcdir, '/untar2.scr ', runsamps$raw1[1], ' ', runsamps$raw2[1], ' ', rundir))
        }
        data$untar[data$run_name==i]='submitted'
    }



    
    
    ##submit calling jobs if untar is done
    numdirs=length(list.dirs(rundir_raw, recursive=FALSE))-1
    numcalled=length(list.files(rundir_done, recursive=FALSE))-1
    
    if (file.exists(untar_done) && (numcalled!=numdirs && runsamps$call[1]!='submitted')) {
        system(paste0('mkdir -p ', rundir_call))
        system(paste0('mkdir -p ', rundir_done))

        system(paste0('sbatch --array=0-', as.character(numdirs), ' --job-name=', name, ' --output=', rundir_done, '/', name, '.%A_%a.out ', srcdir, 'bc_call.scr ', rundir))

        data$call[data$run_name==i]='submitted'
    }



    

    ##gather the fastq for each sample
    for (samp in 1:dim(runsamps)[1]) {
        name=paste0(runsamps$org[samp],'_', runsamps$isolate_num[samp])

        sampdir=paste0(datadir, name)
        fqdir=paste0(sampdir, '/fastqs')
        sampdir_out=paste0(sampdir,'/batch_logs')

        system(paste0('mkdir -p ', sampdir))
        system(paste0('mkdir -p ', fqdir))
        system(paste0('mkdir -p ', sampdir_out))
        
        ##copy the untar so stuff doesn't get confused
        system(paste0('cp ', untar_done, ' ', sampdir, '/'))

        
        fq=paste0(fqdir, '/', name, '.fq')        
        ##gather all the fqs if calling is done
        if (numdirs==numcalled && (!file.exists(fq) && runsamps$fq[samp]!='submitted')) {
            system(paste0('mkdir -p ', sampdir))
            system(paste0('mkdir -p ', fqdir))
            system(paste0('mkdir -p ', sampdir_out))

            barcode=runsamps$barcode[samp]
            system(paste0('sbatch --output=', sampdir_out, '/fq.out ', srcdir, 'bc_fqs.scr ', rundir, ' ', sampdir_fq, ' ', barcode))

            
            fq_done=paste0(fqdir, '/fq_done.txt')

            data$fq[data$isolate_num==runsamps$isolate_num[samp]]='submitted'
            data$call[data$isolate_num==runsamps$isolate_num[samp]]='done'
        }
    }
}

gs_edit_cells(status, input=data, anchor=paste0('A1'))
