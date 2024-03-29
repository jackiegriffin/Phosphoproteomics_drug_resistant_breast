################################################################################################################
## Filename: ssgsea-gui.r
## Created: September 09, 2017
## Author(s): Karsten Krug
##
## Purpose: 
##      - Wrapper to ssGSEA script to perform single sample Gene Set Enrichment analysis.
##          
## Instructions:  
##      - Source the script into a running R-session:
##          - RStudio: open file and press 'Source' in the upper right part of the editor window 
##          - R-GUI: drag and drop this file into an R-GUI window
##      - In order to specify your input files and databases the script will invoke two 
##        Windows file dialogs.
##      - The first dialog lets you choose a folder containing input files in GTC v1.2 format. 
##        The script will loop over all gct files in this directory and run ssGSEA on each file 
##        separately.
##      - The second dialog window lets the user choose a gene set database such as MSigDB. 
##        Some default database can be found in the 'db' subfolder.
##      - MAC users: XQuartz is required for invoking file dialogs
##         
################################################################################################################
rm(list=ls())
script.dir <- dirname(sys.frame(1)$ofile) ## get folder the script is located in
os <- Sys.info()['sysname'] ## determine operating system
if (!require("pacman")) install.packages ("pacman")

## ##########################################################
##  define parameters below:
## ##########################################################

## ssGSEA / PTM-SEA parameters
sample.norm.type    = "rank"              ## "rank", "log", "log.rank", "none" 
weight              = 0.75                ## value between 0 (no weighting) and 1 (actual data counts)
statistic           = "area.under.RES"    ## "Kolmogorov-Smirnov"
output.score.type   = "NES"               ## 'ES' or 'NES'
nperm               = 1e3                ## No. of permutations
min.overlap         = 10                  ## minimal overlap between gene set and data
correl.type         = "z.score"           ## 'rank', 'z.score', 'symm.rank'
par                 = T                   ## use 'doParallel' package?
spare.cores         = 1                   ## No. of cores to leave idle
export.signat.gct   = T                   ## if TRUE gene set GCT files will be exported 
extended.output     = T                   ## if TRUE the GCT files will contain stats on gene set overlaps etc.   

## #####################################################################
##   end paramaters
## - in a perfect world users don't have to worry about the stuff below...
## #####################################################################

## #################################
## directory with gct files
gct.dir.ok=F
while(!gct.dir.ok){
  if(os == 'Windows')
    gct.dir <- choose.dir(default=getwd(), caption = 'Choose directory containing GCT files.')
  else {
    p_load(tcltk)
    gct.dir <- tclvalue(tkchooseDirectory())
  }
  
  if(length(grep('\\.gct$', dir(gct.dir))) > 0)
    gct.dir.ok=T
  
  if(is.na(gct.dir))
    stop('No folder specified! Aborting.')
  }

## directory to write output
out.dir <- gct.dir

## MSigDB
db.ok=F
while(!db.ok){
  
  if(os == 'Windows')
    gene.set.databases = choose.files(default = paste( script.dir, 'db/msigdb/c2.cp.v7.0.symbols.gmt', sep='/' ), caption='Choose gene set database in gmt format. See Broad\'s MSigDB website for details.')
  else
    gene.set.databases = file.choose()
  
  if(length(grep('\\.gmt$', gene.set.databases)) > 0 )
    db.ok=T
  if(length(gene.set.databases)==0)
    stop('No gene set database specified! Aborting.')
  }
## ######################################################################
##                          START
## ######################################################################
#source(paste(script.dir, 'src/ssGSEA2.0.R', sep='/'))
source(file.path(script.dir, 'src/ssGSEA2.0.R'))

## #############################################
## prepare output folder
setwd(out.dir )

date.str <- paste(sub(' .*', '', Sys.time()), sep='_')
dir.create(date.str)
setwd(date.str)


## #############################################
## import signature database
signat.all <- unlist(lapply(gene.set.databases, readLines))
signat.all <- strsplit(signat.all, '\t')
names(signat.all) <- sapply(signat.all, function(x)x[1])
signat.all <- lapply(signat.all, function(x) x[-c(1,2)])

## save parameters used for ssGSEA
param.str = c(
    paste('##', Sys.time()),
    paste('gct.directory:', gct.dir, sep='\t'),
    paste('output.directory:', out.dir, sep='\t'),
    paste('gene.set.database:',gene.set.databases, sep='\t'),
    paste('sample.norm.type:', sample.norm.type, sep='\t'),
    paste('weight:', weight, sep='\t'),
    paste('statistic:', statistic, sep='\t'),
    paste('output.score.type', output.score.type, sep='\t'),
    paste('nperm:', nperm, sep='\t'),
    paste('min.overlap:', min.overlap, sep='\t'),
    paste('correl.type:', correl.type, sep='\t'),
    paste('run.parallel:', par, sep='\t')
   )
writeLines(param.str, con='parameters.txt')


## identify all gct files
gct.files <- dir(gct.dir, pattern='\\.gct$', full.names=T)
names(gct.files) <- paste(  sub('\\.gct$', '', sub('.*/','', gct.files)), 'ssGSEA', sep='_' )

#debug(ssGSEA2)

## #####################################
## loop over gct files and run ssGSEA
for(i in names(gct.files)){


    ## create sub folders if more than one gct file was found
    if(length(gct.files) > 1){
        subdir=sub(',|\\.|:|;|/', '_', i)
        dir.create(subdir)
        setwd(subdir)
    }

    ## ########################################
    ## ssGSEA

    ## input data set
    input.ds <- gct.files[i]

    cat('Running ssSGEA on:', sub('.*/', '', input.ds), '\n\n')

    ## run ssGSEA
    gsea.res <- ssGSEA2(input.ds, gene.set.databases=gene.set.databases, sample.norm.type=sample.norm.type, weight=weight,statistic=statistic, output.score.type = output.score.type, nperm  = nperm, min.overlap  = min.overlap, correl.type = correl.type, output.prefix = paste(i), par=par, 
                        spare.cores=spare.cores, param.file=F, export.signat.gct = export.signat.gct, extended.output = extended.output )

    ## save object
    save(gsea.res, file=paste(i, '.RData', sep=''))

    ## #########################################################
    ##              rank plots
    ## #########################################################
    ## flag to indicate presence of duplicated ids in GCT file
    ## e.g. in case of gene-centric-redundant signature analysis
    dups=F
    if(file.exists( sub('\\.gct', '_unique.gct', input.ds)))
      dups <- T
      
    ## input dataset
    if(dups)
       input.ds <-sub('\\.gct', '_unique.gct', input.ds)
    input.gct <- parse.gctx(input.ds) 
    
    ## gene/site ids
    gn.input <- input.gct@rid
    if(dups)
      gn.input <-  sub('_[0-9]{1,4}$', '', gn.input)
    
    ## sample names
    all.samp <- input.gct@cid
    
    ## expression data only
    input <- input.gct@mat
   
    ## import enrichment scores and p-values
    gsea.score.gct <- parse.gctx(dir('.', pattern=paste( i, '-scores(_[0-9]*x[0-9*]|)', '.gct', sep=''))) 
    gsea.score <- gsea.score.gct@mat
    gsea.pval.gct <- parse.gctx(dir('.', pattern=paste( i,  '-fdr-pvalues(_[0-9]*x[0-9*]|)', sep=''))) 
    gsea.pval <- gsea.pval.gct@mat
    
    ## gene set names
    all.gs <- rownames(gsea.score)

    ## keep only scored signatures
    signat <- signat.all[all.gs]

    ## create sub-folder
    dir.create('rank-plots')
    
    ## loop over gene sets
    for(gs in 1:length(all.gs)){

        gs.name <- all.gs[gs]

        pdf(paste('rank-plots/', make.names( chopString( gsub('\\:|\\/\\\t', ' ', gs.name), nChar=20, add.dots=F)) ,'_2.pdf', sep=''), 9.5, 9.5)
        par(mfrow=c(3, 3))
        for(samp in 1:length(all.samp)){

            ## extract results
            samp.name <- all.samp[samp]

            ## gsea results
            score <- gsea.score[gs.name, samp.name]
            pval <- gsea.pval[gs.name, samp.name]

            ## extract data
            data.expr <- input[, samp.name ]

            valid.idx <- which( !(is.na( data.expr ) | is.infinite(data.expr)) )

            data.expr <- data.expr[ valid.idx ]
            gn <- gn.input[ valid.idx ]

            ## order
            ord.idx <- order(data.expr, decreasing=T)

            ##gn <- row.names(input)[ord.idx]
            gn <- gn[ ord.idx ]
            data.expr <- data.expr[ ord.idx ]


            plot( data.expr, pch=20, col='darkgrey', lwd=4, type='l', xlab='Rank', ylab='Expression', main=paste(gs.name, samp.name, sep='\n'), ylim=range(data.expr), yaxs='i')
						abline(h=0, lty='dashed', lwd=2, col='grey70')

            ## #########################################################
            ##  ptm signatures?
            if(length(grep(';u$|;d$', signat[[gs.name]], value=T)) > 0){

                ## locations
                gsea.tmp.u <- sub(';u$','',grep(';u$', signat[[gs.name]], value=T))
                loc.u <- na.omit(match(gsea.tmp.u, gn))

                gsea.tmp.d <- sub(';d$','',grep(';d$',  signat[[gs.name]], value=T))
                loc.d <- na.omit(match(gsea.tmp.d, gn))

                if(!is.null(loc.u)){

                    rug(loc.u, col='darkred', side=3, lwd=3, ticksize=0.02)
                    points(loc.u, data.expr[loc.u], col=my.col2rgb('darkred',  150), pch=16, cex=2)
                }
                if(!is.null(loc.d)){
                    rug(loc.d, col='darkblue', side=1, lwd=3, ticksize=0.02)
                    points(loc.d, data.expr[loc.d], col=my.col2rgb('darkblue',  150), pch=16, cex=2)

                }
                ## some info
                legend('bottom', legend=c(paste('No. down-regulated in signature:', length(grep(';d$', signat[[gs.name]]))),
                                          paste('No. found in data set:', length(loc.d))
                                          ), inset=.05, bty='n', text.col='darkblue')

                legend('top', legend=c(paste('No. up-regulated in signature:', length(grep(';u$', signat[[gs.name]]))),
                                       paste('No. found in data set:', length(loc.u))
                                       ), inset=.05, bty='n', text.col='darkred')
            } else {## end if signature

                ## ####################################################
                ## regular gene set
                loc <- which(gn %in% signat[[gs.name]])
                rug(loc, col=my.col2rgb('darkred',  50), side=3, lwd=2, ticksize=0.02)
								points(loc, data.expr[loc], col=my.col2rgb('darkred',  150), pch=16, cex=2)

                ## box plot
                loc.quart <- quantile(loc)
                rug(loc.quart, col='darkblue', side=3, lwd=2, ticksize=0.03)
                rect( loc.quart[2], max(data.expr)-0.04*max(data.expr-min(data.expr)), loc.quart[4], max(data.expr), border='darkblue', lwd=2, col=NA )
                rect( loc.quart[1], max(data.expr)-0.02*max(data.expr-min(data.expr)), loc.quart[2], max(data.expr)-0.02*max(data.expr-min(data.expr)), border='darkblue', lwd=2 )
                rect( loc.quart[4], max(data.expr)-0.02*max(data.expr-min(data.expr)), loc.quart[5], max(data.expr)-0.02*max(data.expr-min(data.expr)), border='darkblue', lwd=2 )

                ## some info
                legend('bottom', legend=c(paste('No. in signature:', length( signat[[gs.name]])),
                                          paste('No. found in data set (non-redund.):', sum(signat[[gs.name]] %in% gn)),
                                          paste('No. found in data set (redundant):', length(loc))
                                       ), inset=.05, bty='n', text.col='darkred')

            }
            legend('right', legend=paste('NES=', round(score, 3), ' (p.adj=', round(pval, 5), ')', sep=''), bty='n', inset=.2, cex=1.5)
	}
	par(mfrow=c(1, 1))
	dev.off()
    } ## end loop over gene sets

    if(length(gct.files) > 1)
        setwd('..')
  
    if(dups)
      file.remove(input.ds)
}






