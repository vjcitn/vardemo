workflow task1 {
  call doVariantWorkflow { }
}

task doVariantWorkflow {
  command {
    R -e "BiocManager::install('variants', version = '3.9', update=TRUE, ask=FALSE); \
		library('variants'); \
		file <- system.file('vcf', 'NA06985_17.vcf.gz', package = 'cgdv17'); \
		genesym <- c('TRPV1', 'TRPV2', 'TRPV3'); \
		geneid <- select(org.Hs.eg.db, keys=genesym, keytype='SYMBOL', \
		         columns='ENTREZID'); \
		txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene; \
		seqlevelsStyle(txdb) = 'NCBI'; \
		txdb <- keepSeqlevels(txdb, '17'); \
		txbygene = transcriptsBy(txdb, 'gene'); \
		gnrng <- unlist(range(txbygene[geneid[['ENTREZID']]]), use.names=FALSE); \
		names(gnrng) <- geneid[['SYMBOL']]; \
		param <- ScanVcfParam(which = gnrng, info = 'DP', geno = c('GT', 'cPd')); \
		vcf <- readVcf(file, 'hg19', param); \
		ans = locateVariants(vcf, txdb, AllVariants()); \
		table(mcols(ans)[['LOCATION']]); \
                write.csv(as.data.frame(ans), 'trpvar.csv');"
  }
  runtime {
    docker: "reshg/leonardo-rstudio:develv2"
    }
}
