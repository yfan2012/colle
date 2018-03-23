library('googlesheets')
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library(doParallel))


cl = makeCluster(6)
registerDoParallel(cl, cores=6)

scratch='/scratch/groups/mschatz1/cpowgs/'
srcdir='~/Code/kusama/scripts/'

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

    
    ##untar raw if raw dir isn't there
    if (!dir.exists(paste0(datadir, '/raw'))) {
        system(paste0('mkdir -p ', rawdir))
        
        ##combine minion tars if necessary
        if (is.na(data$raw2[i])) {
            tar=data$raw1[i]
        } else {
            tar=paste0(scratch,name, '_temp.tar')
            system(paste0('tar --concatenate file=',tar,' ', data$raw1[i], ' ', data$raw2[i]))
        }

        ##submit untar script, update sheet, and wait until it's done. Check for completion file every hour. 
        system(paste0('sbatch --output=', batchlogsdir,'/untar_log.txt ', srcdir, 'untar.scr', ' ', tar, ' ', datadir))
        data=gs_read(status)
        data$status[i]='unzipping'
        gs_edit_cells(status, input=data, anchor=paste0('A1'))

        
        wait=0
        while (!file.exists(paste0(datadir, '/untar_done.txt'))){
            Sys.sleep(3600)
            wait=wait+1
            ##If it's taking a weirdly long time, give up on this sample and complain.
            if (wait > 10) {
                data=gs_read(status)
                data$status[i]='ERROR WITH UNZIP'
                gs_edit_cells(status, input=data, anchor=paste0('A1'))
                stop()
            }
        }

        ##update the status now that stuff is unzipped
        data=gs_read(status)
        data$status[i]='unzipped'
        gs_edit_cells(status, input=data, anchor=paste0('A1'))

        if (!is.na(data$raw2[i])) {
            system(paste0('sbatch --output=', batchlogsdir,'/remove_log.txt ' , srcdir, 'remove.scr', ' ', scratch, name, '_temp.tar'))
        }
    }

    
    ##call raw if called isn't there

    if (!dir.exists(calldir)){
        system(paste0('mkdir -p ', calldir))
        system(paste0('mkdir -p ', calllogsdir))
        system(paste0('mkdir -p ', calldonedir))

        
        ##check how large the calling array needs to be, and submit
        numdirs=length(list.dirs(rawdir, recursive=FALSE))-1
        system(paste0('sbatch --array=0-', as.character(numdirs), ' --job-name=', name, ' --output=', calllogsdir, '/', name, '.%A_%a.out ', srcdir , 'call.scr ', datadir))

        ##update to say calling
        data=gs_read(status)
        data$status[i]='calling'
        gs_edit_cells(status, input=data, anchor=paste0('A1'))

        ##wait until all are done
        wait=0
        numdone=length(list.files(calldonedir))
        while (numdone != numdirs+1) {
            Sys.sleep(300)
            wait=wait+1
            numdone=length(list.files(calldonedir))
            if (wait>36) {
                data=gs_read(status)
                data$status[i]='ERROR WITH CALLING'
                gs_edit_cells(status, input=data, anchor=paste0('A1'))
                stop()
            }
        }

        ##update status
        if (numdone != numdirs+1) {
            data=gs_read(status)
            data$status[i]='called'
            data$call[i]=calldir
            gs_edit_cells(status, input=data, anchor=paste0('A1'))
        }
    }



    ##gather fastqs if they aren't there

    if (!dir.exists(fqdir)){
        system(paste0('mkdir -p ', fqdir))

        ##check number of dirs and loop through, looking for the fastq
        numdirs=length(list.dirs(rawdir, recursive=FALSE))-1
        system(paste0('sbatch --output=', batchlogsdir, '/fqs.txt ', srcdir, 'fqs.scr ', datadir, ' ', as.character(numdirs)))

        ##wait until fqs are gathered and placed
        wait=0
        while (!file.exists(paste0(batchlogsdir, '/fqs_done.txt'))) {
            Sys.sleep(300)
            wait=wait+1
            numdone=length(list.files(calldonedir))
            if (wait>36) {
                data=gs_read(status)
                data$status[i]='FQS TAKING A LONG TIME'
                gs_edit_cells(status, input=data, anchor=paste0('A1'))
                stop()
            }
        }
        
        
        ##update fastq location
        if (file.exists(paste0(batchlogsdir, '/fqs_done.txt'))){
            data=gs_read(status)
            data$fq[i]=paste0(fqdir, '/', name, '.fq')
            gs_edit_cells(status, input=data, anchor=paste0('A1'))
        }

    }
    

    
    ##assemble

    if (!dir.exists(canudir)) {
        system(paste0('mkdir -p ', canudir))

        ##start assembly
        system(paste0('bash ', srcdir, 'assemble_bacteria.sh ', datadir))

        ##update to say assembling
        data=gs_read(status)
        data$status[i]='assembling'
        gs_edit_cells(status, input=data, anchor=paste0('A1'))


        ##wait for the assembly to finish. Check every 12 hrs
        wait=-1
        assembly=paste0(datadir, '/canu_assembly/', name, '.contigs.fasta')
        while (!file.exists(assembly)) {
            if (wait > 20) {
                data=gs_read(status)
                data$status[i]='ASSEMBLY ERROR'
                gs_edit_cells(status, input=data[i,], anchor=paste0('A', as.character(i)))
                stop()
            }
            wait=wait+1
            Sys.sleep(43200)
        }
        
        ##update sheet 
        if (!file.exists(assembly)) {
            data=gs_read(status)
            data$status[i]='assembled'
            data$assembly[i]=assembly
            gs_edit_cells(status, input=data[i,], anchor=paste0('A', as.character(i)))
        }
    }


    
    ##polish
    if (!dir.exists(polishdir)) {
        system(paste0('mkidr -p ', polishdir))
        system(paste0('mkdir -p ', polishdir, '/polish_logs'))
        
        
}

    
