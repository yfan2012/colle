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
data=gs_read(status)



##make sure everything is called 
allruns=unique(data$run_code)

foreach(x=1:length(allruns)) %dopar% {
    i=allruns[x]
    
    runsamps=data[data$run_code==i,]
    rundir=paste0(datadir, runsamps$run_name[1])
    rundir_out=paste0(rundir, '/batch_logs')
    rundir_raw=paste0(rundir, '/raw')
    rundir_call=paste0(rundir, '/called')
    rundir_done=paste0(rundir, '/call_done' )
    system(paste0('mkdir -p ', rundir))
    system(paste0('mkdir -p ', rundir_out))
    system(paste0('mkdir -p ', rundir_raw))
    system(paste0('mkdir -p ', rundir_call))
    system(paste0('mkdir -p ', rundir_done))
    
    ##Make sure everything is untared
    if (!identical(runsamps$untar[1],'untared') && !identical(runsamps$untar[1]!='ERROR')) {
        tar1=paste0(epidir, runsamps$tar1[1])
        ##if there ar two tarballs
        if (runsamps$tar2[1]!=NA){
            tar2=paste0(epidir, runsamps$tar2[1])
            system(paste0('sbatch --output=',rundir_out, '/untar.log ', srcdir,'/untar2.scr ', tar1, ' ', tar2, ' ', rundir)) 
        } else {
            system(paste0('sbatch --output=',rundir_out, '/untar.log ', srcdir,'/untar.scr ', tar1, ' ', rundir))
        }

        untar_done=paste0(rundir, '/untar_done.txt')
        
        wait=0
        while (!file.exists(untar_done)) {
            Sys.sleep(3600)
            wait=wait+1
            if (wait>10) {
                data=gs_read(status)
                data$untar[data$run_code==i]='ERROR'
                gs_edit_cells(status, input=data, anchor=paste0('A1'))
                stop()
            }
        }

        
        if (file.exists(untar_done)) {
        data=gs_read(status)
        data$untar[i]='untared'
        gs_edit_cells(status, input=data, anchor=paste0('A1'))
        }
    }


    ##Make sure everything is called
    if (runsamps$untar[1]=='untared' && (!identical(runsamps$call[1],'called') && !identical(runsamps$call[1]!='ERROR'))){
        
        numdirs=length(list.dirs(rundir_raw, recursive=FALSE))-1
        system(paste0('sbatch --array=0-', as.character(numdirs), ' --job-name=', runsamps$run_name[1], ' --output=', rundir_out, '/', name, '.%A_%a.out ', srcdir , 'bc_call.scr ', rundir))
        
        ##wait until all are done
        wait=0
        numdone=length(list.files(rundir_out))
        while (numdone != numdirs+1) {
            Sys.sleep(300)
            wait=wait+1
            numdone=length(list.files(rundir_out))
            if (wait>36) {
                data=gs_read(status)
                data$call[data$run==i]='ERROR'
                gs_edit_cells(status, input=data, anchor=paste0('A1'))
                stop()
            }
        }
        
        ##update status
        if (numdone != numdirs+1) {
            data=gs_read(status)
            data$call[data$run==i]='called'
            gs_edit_cells(status, input=data, anchor=paste0('A1'))
        }
    }


    
    ##gather all the fqs in the proper place
    if (runsamps$call[1]=='called' && (!identical(runsamps$fq[1],'gathered') && !identical(runsamps$fq[1],'ERROR'))){

        ##go through each sample in the run
        foreach(samp=1:dim(runsamps)[1]) %dopar% {
            sampdir=paste0(datadir, runsamps$org[samp], '_', runsamps$isolate[samp])
            sampdir_fq=paste0(sampdir, '/fastqs')
            sampdir_out=paste0(sampdir,'/batch_logs')
            system(paste0('mkdir -p ', sampdir))
            system(paste0('mkdir -p ', sampdir_fq))
            system(paste0('mkdir -p ', sampdir_out))

            barcode=runsamps$barcode[samp]
            print('hi')
            system(paste0('sbatch --output=', sampdir_out, '/fq.out ', srcdir, 'bc_fqs.scr ', rundir, ' ', sampdir_fq, ' ', barcode))

            fq_done=paste0(sampdir_fq, '/fq_done.txt')
            print(fq_done)
            wait=0
            if(FALSE)
            while (!file.exists(fq_done)) {
                Sys.sleep(60)
                wait=wait+1
                if (wait>36) {
                    data=gs_read(status)
                    data$fq[data$isolate==runsamps$isolate[samp]]='ERROR'
                    gs_edit_cells(status, input=data, anchor=paste0('A1'))
                    stop()
                }
            }
            }
            if (file.exists(fq_done)) {
                data=gs_read(status)
                data$fq[data$isolate==runsamps$isolate[samp]]='gathered'
                gs_edit_cells(status, input=data, anchor=paste0('A1'))
            }
        }
    }
}



