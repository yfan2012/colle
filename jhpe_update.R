update_gsheet <- function(data, scratch){
    for (i in 1:dim(data)[1]){
        name=paste0(data$org[i], '_' , data$isolate_num[i])
        rundir=paste0(scratch, data$run_name[i])
        datadir=paste0(scratch, name)
        rawdir=paste0(rundir,'/raw')
        calldir=paste0(rundir,'/called')
        calllogsdir=paste0(rundir,'/call_logs')
        calldonedir=paste0(rundir,'/call_done')
        runlogsdir=paste0(rundir,'/batch_logs')
        batchlogsdir=paste0(datadir, '/batch_logs')
        fqdir=paste0(datadir,'/fastqs')
        canudir=paste0(datadir, '/canu_assembly')
        bamdir=paste0(datadir,'/bams')
        polishdir=paste0(datadir, '/polish')
        mpolishdir=paste0(datadir, '/mpolish')
        
        ##check if untar is done
        if (file.exists(paste0(rundir, '/untar_done.txt'))) {
            data$untar[i]='done'
        }else if (!identical(data$untar[i],'submitted')){
            data$untar[i]='NA'
        }
        
        ##check if calling isdone
        numdirs=length(list.dirs(rawdir, recursive=FALSE))-1
        numcalled=length(list.files(calldonedir, recursive=FALSE))-1
        if (numdirs==numcalled && numdirs!=-1) {
            data$call[i]='done'
        }else if (!identical(data$call[i],'submitted')){
            data$call[i]=NA
        }
        
        ##check if fq exits
        if (file.exists(paste0(fqdir,'/', name, '.fq'))) {
            data$fq[i]='done'
        }else if (!identical(data$fq[i],'submitted')){
            data$fq[i]=NA
        }
        
        ##check if assembly is done
        if (file.exists(paste0(canudir,'/', name, '.contigs.fasta'))) {
            data$assembly[i]='done'
        }else if (!identical(data$assembly[i],'submitted')){
            data$assembly[i]=NA
        }

        ##check if alignment is done
        if (file.exists(paste0(bamdir,'/align_done.txt'))) {
            data$align[i]='done'
        }else if (!identical(data$align[i],'submitted')){
            data$align[i]=NA
        }

        ##check if polish is done
        if (file.exists(paste0(polishdir, '/', name, '.polished.fasta'))) {
            if (file.info(paste0(polishdir, '/', name, '.polished.fasta'))$size!=0) {
                data$polish[i]='done'
            }
        }else if (!identical(data$polish[i], 'submitted')){
            data$polish[i]=NA
        }

        ##check if mpolish is done
        if (file.exists(paste0(mpolishdir, '/', name, '.polished_meth.fasta'))) {
            if (file.info(paste0(mpolishdir, '/', name, '.polished_meth.fasta'))$size!=0) {
                data$mpolish[i]='done'
            }
        }else if (!identical(data$mpolish[i], 'submitted')){
            data$mpolish[i]=NA
        }

    }
    
    return(data)
}


