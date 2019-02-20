workflow task1 {
  call doVariantWorkflow { }
}

task doVariantWorkflow {
  command {
    R -e "BiocManager::install('variants', version = '3.8', update=TRUE, ask=FALSE); \
		library('variants'); \
		file <- system.file('vcf', 'NA06985_17.vcf.gz', package = 'cgdv17'); \
		genesym <- 'ORMDL3'; \
		geneid <- select(org.Hs.eg.db, keys=genesym, keytype='SYMBOL', \
		         columns='ENTREZID'); \
		txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene; \
		seqlevelsStyle(txdb) = 'NCBI'; \
		txdb = keepStandardChromosomes(txdb); \
		txdb <- keepSeqlevels(txdb, '17'); \
		txbygene = transcriptsBy(txdb, 'gene'); \
		gnrng <- unlist(range(txbygene[geneid[['ENTREZID']]]), use.names=FALSE); \
		names(gnrng) <- geneid[['SYMBOL']]; \
		param <- ScanVcfParam(which = gnrng+20000, info = 'DP', geno = c('GT', 'cPd')); \
		vcf <- readVcf(file, 'hg19', param); \
		seqlevels(vcf)[25] = 'MT'; \
		ans = locateVariants(vcf, txdb, AllVariants()); \
		table(mcols(ans)[['LOCATION']]); \
		names(ans) = make.names(names(ans), unique=TRUE); \
		ans = as.data.frame(ans); \
		rownames(ans) = make.names(rownames(ans), unique=TRUE); \
                write.csv(ans, 'trpvar.csv');"
  }
  runtime {
    docker: "waldronlab/bioconductor_devel"
    }
}
