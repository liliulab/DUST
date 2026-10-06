#' Find hypo- and hyper-mutated samples
#' @param input.file.name: A VCF-formatted file to read SnpEff annotated somatic variants. This input file shall contain these fields: Tumor_Sample_Barcode, Chromosome, Start_Position, dbSNP_RS, Reference_Allele, Tumor_Seq_Allele2, FILTER, One_Consequence, Hugo_Symbol, Gene, Feature, ENSP, HGVSc, HGVSp_Short, Amino_acids, Codons, ENSP, RefSeq, Entrez_Gene_Id.
#' @param output.folder: the folder to write the output files
#' @param output.prefix: prefix used to name the output file as prefix.outlier.txt.
#' @return NULL
#' @examples find.outliers(input.file.name='./examples/TCGA.CHOL.somatic.maf.gz', output.folder='./examples/', output.prefix='TCGA.CHOL');
#' @export
find.outliers <- function(input.file.name, output.folder, output.prefix) {
	cat('finding outliers ...\n', '   ', input.file.name, '...'); flush.console();
	outlier.file.name <- paste(output.folder, '/', output.prefix, '.outlier.txt', sep='');
	if(nchar(output.prefix) == 0) {
		outlier.file.name <- paste(working.dir, 'outlier.txt', sep='');
	}

	mut.all <- read.table(input.file.name, sep='\t', header=T, quote='', stringsAsFactors=F);
	cat(nrow(mut.all), 'lines\n'); flush.console();

	if(is.null(mut.all$FILTER)) mut.all$FILTER = 'PASS' 
	mut.filtered <- mut.all[, c('Tumor_Sample_Barcode', 'Chromosome', 'Start_Position', 'dbSNP_RS', 'Reference_Allele', 'Tumor_Seq_Allele2', 'FILTER', 'One_Consequence', 'Hugo_Symbol', 'Gene', 'Feature', 'ENSP', 'HGVSc', 'HGVSp_Short', 'Amino_acids', 'Codons', 'ENSP', 'RefSeq', 'Entrez_Gene_Id')];
	colnames(mut.filtered) <- c('sample', 'chr', 'pos', 'rsid', 'ref', 'alt', 'filter', 'type', 'symbol', 'ensembl.gene', 'cds.id', 'protein.id', 'cds.comp', 'prot.comp', 'aa', 'codon', 'prot.id', 'refseq', 'gene.id');
	mut.filtered <- mut.filtered[which(mut.filtered$type %in% c('synonymous_variant', 'missense_variant', 'stop_gained', 'frameshift_variant', 'inframe_insertion', 'inframe_deletion')), ];
	agg <- aggregate(chr ~ sample, data=mut.filtered, 'length')
	agg <- agg[order(agg$chr), ]
	# lower.cutoff <- m - 3*d;
	# upper.cutoff <- m + 3*d;	
	x <- log(agg$chr, 2)
	m <- mean(x)
	d <- sd(x)
	q1 <- quantile(x, 0.25);
	q3 <- quantile(x, 0.75);
	iqr <- IQR(x)
	lower.cutoff <- 2^(q1 - iqr*1.5);
	upper.cutoff <- 2^(q3 + iqr*1.5);
	lower.outlier <- which(agg$chr < lower.cutoff)
	upper.outlier <- which(agg$chr > upper.cutoff)
	cat('     mean (', round(2^m), '), std (', round(2^d), '), cutoffs (', lower.cutoff, ' - ', upper.cutoff, ')\n'); flush.console();
	cat('     hypo- (', length(lower.outlier), ') ... hyper (', length(upper.outlier), ')\n'); flush.console();
	
	agg$outlier <- 'no';
	agg[lower.outlier, 'outlier'] <- 'hypo';
	agg[upper.outlier, 'outlier'] <- 'hyper';

	write.table(agg, outlier.file.name, sep='\t', quote=F, row.names=F);
}

#' Parse the MAF file to aggregate variants by genes and mutation types
#' @param input.file.name: A VCF-formatted file to read SnpEff annotated somatic variants. This input file shall contain these fields: Tumor_Sample_Barcode, Chromosome, Start_Position, dbSNP_RS, Reference_Allele, Tumor_Seq_Allele2, FILTER, One_Consequence, Hugo_Symbol, Gene, Feature, ENSP, HGVSc, HGVSp_Short, Amino_acids, Codons, ENSP, RefSeq, Entrez_Gene_Id.
#' @param output.folder: the folder to write output files
#' @param output.prefix: prefix used to name the output files. Output files include prefix.error.txt (mutations cannot be mapped), prefix.mut.filtered.txt (somatic variants with matching types, prefix.mut.cnt.txt (somatic variants aggregated by sample, gene and mutational type), prefix.mut.summary.txt (somatic variants aggregated by gene and mutational type across all samples), prefix.symbol_2_cds_id.txt (mapping gene symbols to ensembl transcript ids), and prefix.log.txt (log information).
#' @return NULL
#' @examples parse.aggregated.mut(input.file.name='./examples/TCGA.CHOL.somatic.maf.gz', output.folder='./examples/', output.prefix='TCGA.CHOL')
#' @export
parse.aggregated.mut.2024 <- function(input.file.name, output.folder='', output.prefix='', fa_aa=NULL, fa_bp=NULL, recover=T) {	
	cat('parsing aggregated calls ...\n');
	cat('   ', input.file.name, '...'); flush.console();

	working.dir <- paste(output.folder, '/', sep='');
	error.file.name <- paste(working.dir, output.prefix, '.error.txt', sep='');
	mut.filtered.file.name <- paste(working.dir, output.prefix, '.mut.filtered.txt', sep='');

	mut.all <- read.table(input.file.name, sep='\t', header=T, quote='', stringsAsFactors=F);
	cat(nrow(mut.all), 'lines\n'); flush.console();  
	
	if(is.null(mut.all$FILTER)) mut.all$FILTER = 'PASS' 
	mut.filtered <- mut.all[, c('Tumor_Sample_Barcode', 'Chromosome', 'Start_Position', 'dbSNP_RS', 'Reference_Allele', 'Tumor_Seq_Allele2', 'FILTER', 'One_Consequence', 'Hugo_Symbol', 'Gene', 'Feature', 'ENSP', 'HGVSc', 'HGVSp_Short', 'Amino_acids', 'Codons', 'RefSeq', 'Entrez_Gene_Id', 'DOMAINS')];
	colnames(mut.filtered) <- c('sample', 'chr', 'pos', 'rsid', 'ref', 'alt', 'filter', 'type', 'symbol', 'ensembl.gene', 'cds.id', 'prot.id', 'cds.comp', 'prot.comp', 'aa', 'codon', 'refseq', 'gene.id', 'domain');
	mut.filtered = mut.filtered[which(mut.filtered$filter == 'PASS'), ]
	cat('    mutations with PASS filter: ', nrow(mut.filtered), '\n'); flush.console();	 
	
	mut.filtered <- mut.filtered[which(mut.filtered$type %in% c('synonymous_variant', 'missense_variant', 'stop_gained', 'frameshift_variant', 'inframe_insertion', 'inframe_deletion')), ];
	cat('    mutations of type syn, missense, nonsense, fs-indel, inf-indel', nrow(mut.filtered), '\n'); flush.console();	 ## 136441
	mut.filtered$type <- ifelse(mut.filtered$type == 'synonymous_variant', 'syn', ifelse(mut.filtered$type == 'missense_variant', 'missense', ifelse(mut.filtered$type == 'stop_gained', 'nonsense', ifelse(mut.filtered$type =='frameshift_variant', 'fs', ifelse(mut.filtered$type =='inframe_insertion', 'inf_ins', 'inf_del')))));
	fs.idx <- which(mut.filtered$type == 'fs');
	fs.idx.ins <- intersect(fs.idx, grep('ins', mut.filtered$cds.comp));
	fs.idx.del <- intersect(fs.idx, grep('del', mut.filtered$cds.comp));
	fs.idx.dup <- intersect(fs.idx, grep('dup', mut.filtered$cds.comp));
	if(sum(length(fs.idx.ins), length(fs.idx.del), length(fs.idx.dup)) < length(fs.idx)) {
		cat('!!!!!!!!!!!!!!!!!!! Error !!!!!!!!!!!!! CHECK FS INDELS!!!!!!!!!\n');  flush.console();	
	}
	mut.filtered[fs.idx.ins, 'type'] <- 'fs_ins'
	mut.filtered[fs.idx.del, 'type'] <- 'fs_del'
	mut.filtered[fs.idx.dup, 'type'] <- 'fs_ins'	
	cat('    coding mutations', nrow(mut.filtered ), '; ');  flush.console();
	
	pattern.cds.1 <- stringr::str_match(mut.filtered$cds.comp, 'c\\.([0-9]+)([A-Z]+)>([A-Z]+)')
	pattern.cds.2 <- stringr::str_match(mut.filtered$cds.comp, 'c\\.([0-9]+).*([ins]{3}|[del]{3}|[dup]{3})(.*)')
	pattern.cds.3 <- stringr::str_match(mut.filtered$cds.comp, 'c\\.([0-9]+)_([0-9]+).*([ins]{3}|[del]{3}|[dup]{3})(.*)')
	pattern.prot.1 <- stringr::str_match(mut.filtered$prot.comp, 'p\\.([A-Z*])([0-9]+)([A-Z=*])$')
	pattern.prot.2 <- stringr::str_match(mut.filtered$prot.comp, 'p\\.([A-Z*])([0-9]+)([A-Z])fs(.+)')
	pattern.prot.3 <- stringr::str_match(mut.filtered$prot.comp, 'p\\.([A-Z*])([0-9]+)(.*)')
	pattern.codon <- stringr::str_match(mut.filtered$codon, '([A-Za-z\\-]+)/([A-Za-z\\-]+)')
	pattern.aa <- stringr::str_match(mut.filtered$aa, '([A-Z*\\-]+)/([A-Z*\\-]+)')
	idx = which(is.na(pattern.prot.2[, 2])); pattern.prot.2[idx, c(2, 4)] = pattern.aa[idx, 2:3]
	idx = which(is.na(pattern.prot.3[, 2])); pattern.prot.3[idx, c(2, 4)] = pattern.aa[idx, 2:3]
	mut.filtered$cds.pos <- pattern.cds.1[, 2]
	mut.filtered$cds.ref <- pattern.cds.1[, 3]
	mut.filtered$cds.alt <- pattern.cds.1[, 4]
	idx <- which(is.na(mut.filtered$cds.pos));
	mut.filtered[idx, c('cds.pos', 'cds.ref', 'cds.alt')] <- cbind(pattern.cds.2[idx, 2], mut.filtered[idx, c('ref', 'alt')])
	idx <- which(is.na(mut.filtered$cds.pos));
	mut.filtered[idx, c('cds.pos', 'cds.ref', 'cds.alt')] <- cbind(pattern.cds.3[idx, 2], mut.filtered[idx, c('ref', 'alt')])
	mut.filtered$prot.pos <- pattern.prot.1[, 3]
	mut.filtered$prot.ref <- pattern.prot.1[, 2]
	mut.filtered$prot.alt <- pattern.prot.1[, 4]
	mut.filtered$prot.alt <- ifelse(mut.filtered$prot.alt == '=', mut.filtered$prot.ref, mut.filtered$prot.alt)
	idx <- which(mut.filtered$type %in% c('fs_ins', 'fs_del'));
	mut.filtered[idx, c('prot.pos', 'prot.ref', 'prot.alt')] <- pattern.prot.2[idx, c(3, 2, 4)]
	idx <- which(mut.filtered$type %in% c('inf_ins', 'inf_del'));
	mut.filtered[idx, c('prot.pos', 'prot.ref', 'prot.alt')] <- pattern.prot.3[idx, c(3, 2, 4)]	
	idx <- which(is.na(mut.filtered$prot.pos))
	mut.filtered[idx, c('prot.pos', 'prot.ref', 'prot.alt')] <- cbind(pattern.prot.3[idx, 3], pattern.aa[idx, 2:3])
	mut.filtered$prot.pos = as.numeric(mut.filtered$prot.pos)
	idx <- which(is.na(mut.filtered$prot.pos))
	mut.filtered[idx, 'prot.pos'] <- pattern.prot.1[idx, 3]
	mut.filtered$cds.pos = as.numeric(mut.filtered$cds.pos)
	mut.filtered$codon.ref = toupper(pattern.codon[, 2])
	mut.filtered$codon.alt = toupper(pattern.codon[, 3])	
	
	##run with New Fasta from ENS(Recover BRAF)
	if(!is.null(fa_aa) & !is.null(fa_bp)){	  
	  e.new <- t(apply(mut.filtered[, c('cds.id', 'cds.pos')], 1, function(x) { i<-which(fa_bp$cds.id.2==x[1]); r <- c(); if(length(i)==0) { r <- c('', '', ''); } else { p<-as.numeric(x[2]); s=fa_bp[i, 'cds.seq']; prev=substr(s, p-1, p-1); current=substr(s, p, p); post=substr(s, p+1, p+1); r<-c(prev, current, post);} }))
	  f.new <- unlist(apply(mut.filtered[, c('prot.id', 'prot.pos')], 1, function(x) { i<-which(fa_aa$pep.id.2==x[1]); r <- c(); if(length(i)==0) { r <- c(''); } else { p<-as.numeric(x[2]); s=fa_aa[i, 'pep.seq']; current=substr(s, p, p); r<-c(current); } }))
	  mut.filtered$seq.prev <- e.new[, 1];
	  mut.filtered$seq.current <- e.new[, 2];
	  mut.filtered$seq.post <- e.new[, 3];
	  mut.filtered$prot.current <- f.new;
	  index.error.1 <- which(!mut.filtered$type %in% c('fs_ins', 'fs_del', 'inf_ins', 'inf_del') & mut.filtered$seq.current != mut.filtered$cds.ref);
	  index.error.2 <- which(!mut.filtered$type %in% c('fs_ins', 'fs_del', 'inf_ins', 'inf_del') & mut.filtered$prot.current != mut.filtered$prot.ref);
	  index.error <- unique(c(index.error.1, index.error.2));
	  cat('    New error ', length(index.error), '; ');  flush.console();	  
	}
	
	errors <- c();
	if(length(index.error) > 0) {
	  errors <- mut.filtered[index.error, ]	  
	  mut.filtered <- mut.filtered[-index.error, ]
	}
	
	if(recover){
		##run with old GUST
		e <- t(apply(errors[, c('cds.id', 'cds.pos')], 1, function(x) { i<-which(seq_bp$cds.id.2==x[1]); r <- c(); if(length(i)==0) { r <- c('', '', ''); } else { p<-as.numeric(x[2]); s=seq_bp[i, 'cds.seq']; prev=substr(s, p-1, p-1); current=substr(s, p, p); post=substr(s, p+1, p+1); r<-c(prev, current, post);} }))
		f <- unlist(apply(errors[, c('prot.id', 'prot.pos')], 1, function(x) { i<-which(seq_aa$pep.id.2==x[1]); r <- c(); if(length(i)==0) { r <- c(''); } else { p<-as.numeric(x[2]); s=seq_aa[i, 'pep.seq']; current=substr(s, p, p); r<-c(current); } }))
		errors$seq.prev <- e[, 1];
		errors$seq.current <- e[, 2];
		errors$seq.post <- e[, 3];
		errors$prot.current <- f;
		index.error.1 <- which(!errors$type %in% c('fs_ins', 'fs_del', 'inf_ins', 'inf_del') & errors$seq.current != errors$cds.ref);
		index.error.2 <- which(!errors$type %in% c('fs_ins', 'fs_del', 'inf_ins', 'inf_del') & errors$prot.current != errors$prot.ref);
		index.error <- unique(c(index.error.1, index.error.2));
		cat('    Old Gust error ', length(index.error), '; ');  flush.console();
		
		if(length(index.error) > 0) {
		  errors <- errors[index.error, ]		  
		  mut.filtered <- rbind(mut.filtered,errors[-index.error, ])
		}else{
		  mut.filtered <- rbind(mut.filtered,errors)
		  errors<-c()
		}	
	}
	
	if (length(errors) > 0) {
		agg.error <- aggregate(chr ~ symbol + pos + ensembl.gene + cds.id + prot.id, data=errors, 'length');
		agg.error <- agg.error[order(-agg.error$chr), ]
		agg.error <- aggregate(chr ~ symbol +  ensembl.gene + cds.id + prot.id, data=errors, 'length');
		agg.error <- agg.error[order(-agg.error$chr), ]	  
		cat('    unmatched mutations', length(index.error), 'in', length(unique(agg.error$symbol)), 'genes. max cnt in a gene :', agg.error[1, 'symbol'], ':', agg.error[1, 'chr']); flush.console();
	}
	write.table(errors, error.file.name, sep='\t', row.names=F, quote=F);
  
	cat(' with matching seq', nrow(mut.filtered), '\n'); flush.console();	## 66649 ##328697
	if(nrow(mut.filtered) > 0) {
		mut.filtered$mut_type <- paste(mut.filtered$cds.ref, mut.filtered$cds.alt, sep='-');
		mut.filtered$mut_type <- ifelse(mut.filtered$mut_type == 'C-T' & mut.filtered$seq.post == 'G', 'C-T-h', ifelse(mut.filtered$mut_type == 'G-A' & mut.filtered$seq.prev == 'C', 'G-A-h', mut.filtered$mut_type));
		mut.filtered$mut_type <- ifelse(mut.filtered$type == 'fs_ins', 'fs_ins', mut.filtered$mut_type);
		mut.filtered$mut_type <- ifelse(mut.filtered$type == 'fs_del', 'fs_del', mut.filtered$mut_type);
		mut.filtered$mut_type <- ifelse(mut.filtered$type == 'inf_ins', 'inf_ins', mut.filtered$mut_type);
		mut.filtered$mut_type <- ifelse(mut.filtered$type == 'inf_del', 'inf_del', mut.filtered$mut_type);
		mut.filtered <- merge(mut.filtered, annotation, by.x='ensembl.gene', by.y='ensembl.gene.id', all.x=T, all.y=F);
		mut.filtered$gene.id <- as.integer(as.character(mut.filtered$gene.id.y));
		mut.filtered$symbol <- mut.filtered$symbol.y;	
		write.table(mut.filtered, mut.filtered.file.name, sep='\t', row.names=F, quote=F);
	}
}

#' Annotate mutations by CDD region. It requires the output.folder contains a file from running the parse.aggregated.mut() function. That is, prefix.mut.filtered.txt file must exist.
#' @param output.folder: the folder to write output files
#' @param output.prefix: prefix used to name the output file as prefix.mut.region.txt. 
#' @return NULL
#' @examples annotate.by.cdd.region(output.folder='./examples/', output.prefix='TCGA.CHOL')
#' @export
annotate.by.cdd.region = function(output.folder='', output.prefix='') {
	cat('annotating mutations by CDD regions ...\n'); flush.console();
	working.dir <- paste(output.folder, '/', sep='');
	cat('    folder', working.dir, '...');
	mut.file.name <- paste(working.dir, output.prefix, '.mut.filtered.txt', sep='');
	mut.region.file.name <- paste(working.dir, output.prefix, '.mut.region.txt', sep='');
	mut.filtered = read.table(mut.file.name, sep='\t', header=T, stringsAsFactors=F)
	cat('   ', nrow(mut.filtered), 'mutations \n')

	mm.prot.idx = c()
	mm.rna.idx = c()
	mut.cdd = c()
	for(i in 1:nrow(mut.filtered)) {   ## i = which(mut.filtered$symbol == 'BRAF' & nchar(mut.filtered$domain) > 5)[1]
		if(i %% 1000 == 0) cat(i, '...\n'); flush.console()
		one = mut.filtered[i, ]
		ensp = one$prot.id
		nm = gsub('\\..+', '', one$refseq)
		prot.pos = one$prot.pos
		type = one$type
		cdd = cdd.region.mapped[which(cdd.region.mapped$ensp == ensp & cdd.region.mapped$nm == nm), ]
		if(nrow(cdd) == 0) cdd = cdd.region.mapped[which(cdd.region.mapped$nm == nm & cdd.region.mapped$refseq.rna == nm), ]
		if(nrow(cdd) == 0) cdd = cdd.region.mapped[which(cdd.region.mapped$ensg == one$ensembl.gene & cdd.region.mapped$refseq.rna == nm), ]
		if(nrow(cdd) == 0) {
			cdd = cdd.region.mapped[which(cdd.region.mapped$ensg == one$ensembl.gene), ]
			temp = gsub('NP_', '', cdd$refseq)
			temp.sel = temp[which.min(as.numeric(temp))]
			cdd = cdd[which(cdd$refseq == paste0('NP_', temp.sel)), ]
		}
		cdd = unique(cdd)
		cdd.by.prot = cdd[which(cdd$start <= prot.pos & cdd$end >= prot.pos), ]
		if(nrow(cdd.by.prot) > 0) {
			aa.coord = prot.pos - cdd.by.prot$start + 1
			aa = substr(cdd.by.prot$seq.aa, aa.coord, aa.coord)
			idx = which(aa == one$prot.ref)
			if(length(idx) != length(aa)) {
				mm.prot.idx = c(mm.prot.idx , i)
			}
			matched = cdd.by.prot[idx, ]
			matched$mut.aa.coord = aa.coord[idx]
			if(nrow(matched) > 0) {
				matched$type = type
				bp.coord = matched$mut.aa.coord * 3
				bp = substr(matched$seq.bp, bp.coord-2, bp.coord)
				idx.2 = which(bp == one$codon.ref)
				if(type %in% c('syn', 'missense', 'nonsense')) {
					if(length(idx.2) != length(bp)) {
						mm.rna.idx = c(mm.rna.idx , i)
					}
					matched.2 = matched[idx.2, ]
					matched.2$mut.codon.coord = bp.coord[idx.2]
				} else {  ## no checking on codon matching for frameshift mut.
					matched.2 = matched
					matched.2$mut.codon.coord = bp.coord
				}

				if(nrow(matched.2) > 0) {
					matched.2$mut.idx = i
					mut.cdd = rbind(mut.cdd, matched.2)
				}
			}
		}		
	}
	cat(' in', nrow(mut.cdd), 'CDD regions');   ## 41487
	cat('     with mismatches found in', length(mm.prot.idx), 'proteins and', length(mm.rna.idx), ' mRNAs\n');   ##  134
	x = substr(mut.cdd$seq.aa, mut.cdd$mut.aa.coord, mut.cdd$mut.aa.coord)
	y = mut.filtered[mut.cdd$mut.idx, 'prot.ref']
	cat('   ', length(which(x != y)), ' mismatched amino acids')  ## 0
	x = substr(mut.cdd$seq.bp, mut.cdd$mut.codon.coord-2, mut.cdd$mut.codon.coord)
	y = mut.filtered[mut.cdd$mut.idx, 'codon.ref']
	cat(' and', length(which(x != y)), 'mismatched nucleotides\n')  ## 715 because of indels
	z = mut.filtered[mut.cdd$mut.idx, 'type']
	table(z[which(x == y)])
	table(z[which(x != y)])
	agg = as.data.frame(table(mut.cdd$mut.idx))
	agg = agg[order(-agg$Freq), ]
	cat(nrow(agg), 'unique mutations in CDD regions\n');  ## 29403 mutations in CDD regions
	# length(unique(grep('Pfam', mut.filtered$domain), grep('SMART', mut.filtered$domain), grep('TIGR', mut.filtered$domain)))  ## 29603
	agg = as.data.frame(table(mut.cdd$db_xref))
	agg = agg[order(-agg$Freq), ]
	cat(nrow(agg), 'unique CDD regions mutated.\n');  ## 8342 CDD regions mutated
	head(agg)
	table(mut.cdd[which(mut.cdd$db_xref %in% agg[1:10, 'Var1']), 'note'])
	
	write.table(mut.cdd, mut.region.file.name, sep='\t', row.names=F, quote=F);
}

#' Calculate the expected frequencies for CDD regions. It requires the output.folder contains a file from running the annotate.by.cdd.region() function. That is, prefix.mut.region.txt file must exist.
#' @param output.folder: the folder to write output files
#' @param output.prefix: prefix used to name the output file as prefix.exp_seq.by_type.region.txt. 
#' @return NULL
#' @examples calculate.expected.by.region(output.folder='./examples/', output.prefix='TCGA.CHOL')
#' @export
calculate.expected.by.region = function(output.folder=output.folder, output.prefix=output.prefix) {
	cat('calculating expected frequencies for CDD regions ...\n'); flush.console();
	working.dir <- paste(output.folder, '/', sep='');
	mut.region.file.name <- paste(working.dir, output.prefix, '.mut.region.txt', sep='');
	region.expected.file.name <- paste(working.dir, output.prefix, '.exp_seq.by_type.region.txt', sep='');

	mut.cdd <- read.table(mut.region.file.name, sep='\t', quote='', header=T, stringsAsFactors=F);
	results.codon_cat$ID <- paste(results.codon_cat$prev, results.codon_cat$codon, results.codon_cat$post, sep='');
	cdd.uniq = unique(mut.cdd$db_xref)
	cat('   ', length(cdd.uniq), 'unique CDD db_xref\n');  ## 8342
	cat('    processing...'); flush.console()
	results.seq_cat <- c();
	for(i in 1:length(cdd.uniq)) {
		if(i %% 1000 == 0) {
			cat(' ', i, '...'); flush.console();
		}
		acc = cdd.uniq[i]
		cdd.this = unique(mut.cdd[which(mut.cdd$db_xref == acc), c('db_xref', 'refseq', 'start', 'end', 'seq.aa', 'seq.bp')])
		exp_syn.list = exp_nonsyn.list = exp_nonsense.list = list()
		for(g in 1:nrow(cdd.this)) {
			seq.bp <- cdd.this[g, 'seq.bp'];
			seq.bp = seqinr::s2c(seq.bp)
			seq.length <- length(seq.bp);
			pos.prev <- seq(3, seq.length, 3);
			pos.post <- seq(4, seq.length, 3);
			codon.cnt <- seq.length/3;
			pos.codon.1 <- seq(1, codon.cnt*3, 3);
			pos.codon.2 <- seq(2, codon.cnt*3, 3);
			pos.codon.3 <- seq(3, codon.cnt*3, 3);
			prev <- c('%', seq.bp[pos.prev]);
			post <- c(seq.bp[pos.post], '%');
			codon <- paste(seq.bp[pos.codon.1], seq.bp[pos.codon.2], seq.bp[pos.codon.3], sep='');
			prev <- prev[1:length(codon)]
			post <- post[1:length(codon)]
			this <- data.frame(prev, codon, post, ID=paste(prev, codon, post, sep=''), stringsAsFactors=F);
			result <- merge(results.codon_cat, this, by='ID', all.x=F, all.y=T);
			result <- result[which(rowSums(result[, 5:7]) > 0), ]
			exp_syn <- aggregate(cnt_syn ~ mut_type, data=result, 'sum'); 
			exp_nonsyn <- aggregate(cnt_nonsyn ~ mut_type, data=result, 'sum'); 
			exp_nonsense <- aggregate(cnt_nonsense ~ mut_type, data=result, 'sum'); 
			exp_syn.list[[g]] = exp_syn
			exp_nonsyn.list[[g]] = exp_nonsyn
			exp_nonsense.list[[g]] = exp_nonsense			
		}
		exp_syn = Reduce(function(x, y) merge(x, y, by='mut_type', all=T), exp_syn.list)
		exp_nonsyn = Reduce(function(x, y) merge(x, y, by='mut_type', all=T), exp_nonsyn.list)
		exp_nonsense = Reduce(function(x, y) merge(x, y, by='mut_type', all=T), exp_nonsense.list)
		if(nrow(cdd.this) > 1) {
			exp_syn = data.frame(mut_type=exp_syn$mut_type, cnt_syn=rowSums(exp_syn[, 2:ncol(exp_syn)]), stringsAsFactors=F)
			exp_nonsyn = data.frame(mut_type=exp_nonsyn$mut_type, cnt_nonsyn=rowSums(exp_nonsyn[, 2:ncol(exp_nonsyn)]), stringsAsFactors=F)
			exp_nonsense = data.frame(mut_type=exp_nonsense$mut_type, cnt_nonsense=rowSums(exp_nonsense[, 2:ncol(exp_nonsense)]), stringsAsFactors=F)
		}
		exp <- merge(merge(exp_syn, exp_nonsyn, by='mut_type'), exp_nonsense, by='mut_type');
		exp <- exp[order(exp$mut_type), ]
		exp$acc <- acc;
		results.seq_cat <- rbind(results.seq_cat, exp);
	}
	cat('\n');
	cat('   ', nrow(results.seq_cat), 'expected sequence by mut_type \n');  ## 116698       5          
	write.table(results.seq_cat, region.expected.file.name, sep='\t', row.names=F, quote=F);
}

#' Organize a collection of variants in a gene into counts 
#' @param x: a data frame with four columns: name, type, mut_type, cnt.
#' @return: a data frame with counts of variants grouped by types
organize.obs <- function(x) {
	colnames(x) <- c('name', 'type', 'mut_type', 'cnt');
	mt <- unique(x$mut_type);
	syn.sum <- data.frame(mut_type=mt, sample=rep(0, length(mt)));
	missense.sum <- syn.sum;
	nonsense.sum <- syn.sum;
	fs.sum <- syn.sum;
	inframe.sum <- syn.sum;
	s <- which(x$type == 'syn');
	m <- which(x$type == 'missense');
	n <- which(x$type == 'nonsense');
	fi <- which(x$type == 'fs_ins');
	fd <- which(x$type == 'fs_del');
	ii <- which(x$type == 'inf_ins');
	id <- which(x$type == 'inf_del');
	if(length(s) > 0) {
		syn.sum <- x[s, c('mut_type', 'cnt')];
	}
	if(length(m) > 0) {
		missense.sum <- x[m, c('mut_type', 'cnt')];
	}
	if(length(n) > 0) {
		nonsense.sum <- x[n, c('mut_type', 'cnt')];
	}
	if(length(fi) + length(fd) > 0) {
		fs.sum <- x[c(fi, fd), c('mut_type', 'cnt')];
	}
	if(length(ii) + length(id) > 0) {
		inframe.sum <- x[c(ii, id), c('mut_type', 'cnt')];
	}
	colnames(syn.sum)[2] <- 'cnt.syn';
	colnames(missense.sum)[2] <- 'cnt.missense';
	colnames(nonsense.sum)[2] <- 'cnt.nonsense';
	colnames(fs.sum)[2] <- 'cnt.fs';
	colnames(inframe.sum)[2] <- 'cnt.inframe';
	sum.obs <- Reduce(function(...) merge(..., all=T, by='mut_type'), list(syn.sum, missense.sum, nonsense.sum, fs.sum, inframe.sum))
	colnames(sum.obs) <- c('mut_type', 'syn.obs', 'missense.obs', 'nonsense.obs', 'fs.obs', 'inframe.obs');
	sum.obs[is.na(sum.obs)] <- 0;
	sum.obs$mut_type <- as.character(sum.obs$mut_type);
	return(sum.obs);	
}

#' Combine reciprocal mutation types
#' @param obs_exp: A data frame with observed and expected counts of variants
#' @return A data frame with observed and expected counts of variants after combining reciprocal mutation types
condense.type <- function(obs_exp) {
	AC.TG <- colSums(obs_exp[which(obs_exp$mut_type %in% c('A-C', 'T-G')), -1], na.rm=T);
	AG.TC <- colSums(obs_exp[which(obs_exp$mut_type %in% c('A-G', 'T-C')), -1], na.rm=T);
	AT.TA <- colSums(obs_exp[which(obs_exp$mut_type %in% c('A-T', 'T-A')), -1], na.rm=T);
	CA.GT <- colSums(obs_exp[which(obs_exp$mut_type %in% c('C-A', 'G-T')), -1], na.rm=T);
	CG.GC <- colSums(obs_exp[which(obs_exp$mut_type %in% c('C-G', 'G-C')), -1], na.rm=T);
	GA.CT <- colSums(obs_exp[which(obs_exp$mut_type %in% c('G-A', 'C-T')), -1], na.rm=T);
	CTh.GAh <- colSums(obs_exp[which(obs_exp$mut_type %in% c('C-T-h', 'G-A-h')), -1], na.rm=T);
	fs_ins <- colSums(obs_exp[which(obs_exp$mut_type %in% c('fs_ins')), -1]);
	fs_del <- colSums(obs_exp[which(obs_exp$mut_type %in% c('fs_del')), -1]);
	inf_ins <- colSums(obs_exp[which(obs_exp$mut_type %in% c('inf_ins')), -1]);
	inf_del <- colSums(obs_exp[which(obs_exp$mut_type %in% c('inf_del')), -1]);
	
	result <- rbind(AC.TG, AG.TC, AT.TA, CA.GT, CG.GC, GA.CT, CTh.GAh, fs_ins, fs_del, inf_ins, inf_del);
	result <- data.frame(mut_type=rownames(result), result, stringsAsFactors=F);
	result[is.na(result)] <- 0;
	return(result);
}

#' optimization function to estimate independent probabilities of a specific mutation type 
prob <- function(theta, obs_exp, compare) {
	mut_type <- unique(obs_exp$mut_type);
	if(compare == 'missense') {
		mut_type <- mut_type[which(mut_type != 'fs')];
	}
	p <- 0;
	for(t in mut_type) {
		obs_syn <- obs_exp[which(obs_exp$mut_type == t), 'syn.obs'];
		exp_syn <- obs_exp[which(obs_exp$mut_type == t), 'syn.exp'];
		obs_compare <- 0;
		exp_compare <- 0;
		obs_total <- 0;
		if(compare == 'missense') {
			obs_compare <- obs_exp[which(obs_exp$mut_type == t), 'missense.obs'];
			exp_compare <- obs_exp[which(obs_exp$mut_type == t), 'missense.exp'];
		} else if(compare == 'nonsense') {
			obs_compare <- sum(obs_exp[which(obs_exp$mut_type == t), c('nonsense.obs', 'fs.obs')]);  ## mut_type=='fs', nonsense.obs=0; otherwise, fs.obs=0;
			exp_compare <- obs_exp[which(obs_exp$mut_type == t), 'nonsense.exp'];
		} else {
			obs_compare <- sum(obs_exp[which(obs_exp$mut_type == t), c('missense.obs', 'nonsense.obs', 'fs.obs', 'inframe.obs')]);
			exp_compare <- sum(obs_exp[which(obs_exp$mut_type == t), c('missense.exp', 'nonsense.exp')]);
		}
		obs_total <- obs_syn + obs_compare;
		if(exp_syn == 0 | exp_compare == 0) {
			pp <- 0;
		} else {
			# pp <- lfactorial(obs_total) - lfactorial(obs_syn) - lfactorial(obs_compare) + obs_syn*log(exp_syn) + obs_compare*log(theta*exp_compare) - obs_total*log(exp_syn + theta * exp_compare);
			pp <- lfactorial(obs_total) - lfactorial(obs_syn) - lfactorial(obs_compare) + obs_syn*log(exp_syn) + obs_compare*theta + obs_compare*log(exp_compare) - obs_total*log(exp_syn + exp(theta)*exp_compare);
			# pp <- lfactorial(obs_total) - lfactorial(obs_inframe) - lfactorial(obs_fs) + obs_inframe*log(exp_inframe) + obs_fs*psi + obs_fs*log(exp_fs) - obs_total*log(exp_inframe + exp(psi)*exp_fs);
		}
		p <- p + pp;
		# cat(t, '\t', pp, '\t', p, '\n');
	}
	return(-p);
}

#' optimization function to estimate joint probabilities of different types of variants. Indels are treated independently from substitutions.
prob.joint <- function(par, obs_exp) {
#	cat(par);
	phi <- par[1];
	psi <- par[2];
	if(phi < -5 | phi > 5 | psi < -5 | psi > 5) {
		return(NA);
	}
	mut_type <- unique(obs_exp$mut_type);
	p <- 0;
	for(t in mut_type) {
#		cat(t, ':' );
		obs_syn <- obs_exp[which(obs_exp$mut_type == t), 'syn.obs'];
		exp_syn <- obs_exp[which(obs_exp$mut_type == t), 'syn.exp'];
		obs_missense <- obs_exp[which(obs_exp$mut_type == t), 'missense.obs'];
		exp_missense <- obs_exp[which(obs_exp$mut_type == t), 'missense.exp'];
		obs_nonsense <- obs_exp[which(obs_exp$mut_type == t), 'nonsense.obs'];
		exp_nonsense <- obs_exp[which(obs_exp$mut_type == t), 'nonsense.exp'];
		exp_syn <- ifelse(exp_syn == 0, 1, exp_syn);
		exp_missense <- ifelse(exp_missense == 0, 1, exp_missense);
		exp_nonsense <- ifelse(exp_nonsense == 0, 1, exp_nonsense);
		
		pp <- 0;
		if(t == 'inframe') {
		} else if(t == 'fs') {	
			obs_fs <- obs_exp[which(obs_exp$mut_type == t), 'fs.obs'];
			obs_inframe <- obs_exp[which(obs_exp$mut_type == 'inframe'), 'inframe.obs'];
			if(obs_fs > 0 | obs_inframe > 0) {
				gene.length <- sum(obs_exp[which(obs_exp$mut_type == t), c('syn.exp', 'missense.exp', 'nonsense.exp')])/3;
				exp_fs <- gene.length*2/3;  ## expected fs is the length of the protein
	#			exp_fs <- obs_exp[which(obs_exp$mut_type == t), 'nonsense.exp']; 
				exp_inframe <- gene.length*1/3;  ## expected fs is the length of the protein
				obs_total <- obs_inframe + obs_fs;
#				pp <- lfactorial(obs_total) - lfactorial(obs_syn) - lfactorial(obs_fs) + obs_syn*log(exp_syn) + obs_fs*log(psi*exp_fs) - obs_total*log(exp_syn + psi*exp_fs);
#				pp <- lfactorial(obs_total) - lfactorial(obs_syn) - lfactorial(obs_fs) + obs_syn*log(exp_syn) + obs_fs*psi + obs_fs*log(exp_fs) - obs_total*log(exp_syn + exp(psi)*exp_fs);
				pp <- lfactorial(obs_total) - lfactorial(obs_inframe) - lfactorial(obs_fs) + obs_inframe*log(exp_inframe) + obs_fs*psi + obs_fs*log(exp_fs) - obs_total*log(exp_inframe + exp(psi)*exp_fs);
			}
#			pp=0;
		} else {
			obs_total <- obs_syn + obs_missense + obs_nonsense;
			if(obs_total > 0) {
#				pp <- lfactorial(obs_total) - lfactorial(obs_syn) - lfactorial(obs_missense) - lfactorial(obs_nonsense) + obs_syn*log(exp_syn) + obs_missense*log(phi*exp_missense) + obs_nonsense*log(psi*exp_nonsense) - obs_total*log(exp_syn + phi*exp_missense + psi*exp_nonsense);
				pp <- lfactorial(obs_total) - lfactorial(obs_syn) - lfactorial(obs_missense) - lfactorial(obs_nonsense) + obs_syn*log(exp_syn) + obs_missense*phi + obs_missense*log(exp_missense) + obs_nonsense*psi + obs_nonsense*log(exp_nonsense) - obs_total*log(exp_syn + exp(phi)*exp_missense + exp(psi)*exp_nonsense);
			}
		}
		p <- p + pp;
#		cat(t, '\t', pp, '\t', p, '\n');
	}
#	cat(' ', -p, '\n');
	return(-p);
}

#' optimization function to estimate joint probabilities of different types of variants including indels.
prob.joint.indel <- function(par, obs_exp) {
	# cat(par);
	phi <- par[1];
	psi <- par[2];
	if(phi < -5 | phi > 5 | psi < -5 | psi > 5) {
		return(NA);
	}
	mut_type <- unique(obs_exp$mut_type);
	mut_type <- mut_type[which(mut_type != 'NA-NA')]
	p <- 0;
	for(t in mut_type) {
		# cat(t, ':' );
		obs_syn <- obs_exp[which(obs_exp$mut_type == t), 'syn.obs'];
		exp_syn <- obs_exp[which(obs_exp$mut_type == t), 'syn.exp'];
		obs_missense <- obs_exp[which(obs_exp$mut_type == t), 'missense.obs'];
		exp_missense <- obs_exp[which(obs_exp$mut_type == t), 'missense.exp'];
		obs_nonsense <- obs_exp[which(obs_exp$mut_type == t), 'nonsense.obs'];
		exp_nonsense <- obs_exp[which(obs_exp$mut_type == t), 'nonsense.exp'];
		exp_syn <- ifelse(exp_syn == 0, 1, exp_syn);
		exp_missense <- ifelse(exp_missense == 0, 1, exp_missense);
		exp_nonsense <- ifelse(exp_nonsense == 0, 1, exp_nonsense);
		
		pp <- 0;
		if(t %in% c('inf_ins', 'inf_del')) {
		} else if(t %in% c('fs_ins', 'fs_del')) {	
			gene.length <- sum(exp_syn, exp_missense, exp_nonsense)/3;
			exp_fs <- gene.length * 2/3;
			exp_inframe <- gene.length * 1/3;
			obs_fs <- obs_exp[which(obs_exp$mut_type == t), 'fs.obs'];
			obs_inframe <- 0; 
			if(t == 'fs_ins') {
				obs_inframe <- obs_exp[which(obs_exp$mut_type == 'inf_ins'), 'inframe.obs'];
				exp_fs <- exp_fs * 0.34;
				exp_inframe <- exp_inframe * 0.34;
			} else {
				obs_inframe <- obs_exp[which(obs_exp$mut_type == 'inf_del'), 'inframe.obs'];
				exp_fs <- exp_fs * 0.66;
				exp_inframe <- exp_inframe * 0.66;
			}
			obs_total <- obs_inframe + obs_fs;
			if(obs_total > 0) {
				pp <- lfactorial(obs_total) - lfactorial(obs_inframe) - lfactorial(obs_fs) + obs_inframe*log(exp_inframe) + obs_fs*psi + obs_fs*log(exp_fs) - obs_total*log(exp_inframe + exp(psi)*exp_fs);
			}
		} else {
			obs_total <- obs_syn + obs_missense + obs_nonsense;
			if(obs_total > 0) {
				pp <- lfactorial(obs_total) - lfactorial(obs_syn) - lfactorial(obs_missense) - lfactorial(obs_nonsense) + obs_syn*log(exp_syn) + obs_missense*phi + obs_missense*log(exp_missense) + obs_nonsense*psi + obs_nonsense*log(exp_nonsense) - obs_total*log(exp_syn + exp(phi)*exp_missense + exp(psi)*exp_nonsense);
			}
		}
		p <- p + pp;
		# cat(t, '\t', pp, '\t', p, '\n');
	}
	# cat(' ', -p, '\n');
	return(-p);
}

#' Estimate selection coefficients of missense mutations and truncating mutations
#' @param obs_exp: A data frame with observed and expected counts
#' @param joint: A boolean flag to estimate joint probabilities
#' @param indel: A boolean flag to consider indel jointly with or independently from substitutions 
#' @return A vector with observed counts, end values and estimated probabilities
compute.selection <- function(obs_exp, joint=FALSE, indel=FALSE) {
	obs_exp[is.na(obs_exp)] <- 0;
	sums <- colSums(obs_exp[-1]);
	sum.obs <- sum(sums[c('syn.obs', 'missense.obs', 'nonsense.obs', 'fs.obs', 'inframe.obs')]);
	sum.obs.disabling <- sum(sums[c('nonsense.obs', 'fs.obs')]);
	sum.exp <- sum(sums[c('syn.exp', 'missense.exp', 'nonsense.exp')]);
	opt.missense <- c(NA, NA);
	opt.nonsense <- c(NA, NA);
	opt.prot <- c(NA, NA);
	if(sum.obs >= 3 & sum.exp > 0) {
		if(joint) {
			opt.joint <- c();
			if(indel) {
				opt.joint <- unlist(optim(par=c(0, 0), fn=prob.joint.indel, obs_exp=obs_exp)); 
			} else {
				opt.joint <- unlist(optim(par=c(0, 0), fn=prob.joint, obs_exp=obs_exp)); 
			}
			if(sum(sums[c('syn.obs', 'missense.obs')]) > 0) opt.missense <- as.numeric(opt.joint[c(1,3)]); 
			if(sum(sums[c('syn.obs', 'nonsense.obs', 'fs.obs')]) > 0) opt.nonsense <- as.numeric(opt.joint[c(2,3)]);
			opt.prot <- unlist(optimize(prob, interval=c(-5, 5), obs_exp, 'prot'));
		} else {
			a <- sums[1]; b <- sums[2]; c <- sums[3]; d <- sums[4]; e <- sums[5];
			if(b > 1 & (a + b) >= 3) {
				opt.missense <- unlist(optimize(prob, interval=c(-5, 5), obs_exp, 'missense')); 
			}
			if((c + d) > 1 & (a + c + d) >= 3) {
				opt.nonsense <- unlist(optimize(prob, interval=c(-5, 5), obs_exp, 'nonsense')); 
			}
			if((b + c + d) > 1 & (a + b + c + d) >= 3) {
				opt.prot <- unlist(optimize(prob, interval=c(-5, 5), obs_exp, 'prot')); 
			}
		}
	}
	opt <- c(sums[1:5], opt.missense, opt.nonsense, opt.prot);
	return(opt);
}

#' Estimate selection coefficients of missense mutations and truncating mutations for CDD regions. It requires the output.folder contains files from running the parse.aggregated.mut(), annotate.by.cdd.region(), and calculate.expected.by.region() functions. That is, prefix.mut.filtered.txt, prefix.mut.region.txt, and prefix.exp_seq.by_type.region.txt files must exist.
#' @param output.folder: the folder to write output files
#' @param output.prefix: prefix used to name the output file as prefix.selection.region.txt. 
#' @param N.sim: the number of simulations to perform to estimate false positives (default: 1000).
#' @return NULL
#' @examples compute.selection.by.region(output.folder='./examples/', output.prefix='TCGA.CHOL')
#' @export
compute.selection.by.region = function(output.folder=output.folder, output.prefix=output.prefix, N.sim=1000) {
	cat('computing selection coefficients for CDD regions ...\n');  flush.console();
	working.dir <- paste(output.folder, '/', sep='');
	mut.region.file.name <- paste(working.dir, output.prefix, '.mut.region.txt', sep='');
	mut.file.name <- paste(working.dir, output.prefix, '.mut.filtered.txt', sep='');
	region.expected.file.name <- paste(working.dir, output.prefix, '.exp_seq.by_type.region.txt', sep='');
	selection.region.file.name <- paste(working.dir, output.prefix, '.selection.region.txt', sep='');

	mut.cdd <- read.table(mut.region.file.name, sep='\t', quote='', header=T, stringsAsFactors=F);
	results.seq_cat <- read.table(region.expected.file.name, sep='\t', header=T, stringsAsFactors=F);
	mut.filtered = read.table(mut.file.name, sep='\t', header=T, stringsAsFactors=F)

	cdd.uniq = unique(mut.cdd$db_xref)
	cat('   ', length(cdd.uniq), 'mutated regions\n');  ## 8342
	result.opt = c()
	result.opt.perm = list()
	cat('    processing...'); flush.console()
	for(i in 1:length(cdd.uniq)) {
		if(i %% 1000 == 0) cat(' \n...', i, '...'); flush.console();
		
		acc = cdd.uniq[i]
		cdd.this = mut.cdd[which(mut.cdd$db_xref == acc), ]
		# dim(cdd.this); table(cdd.this$symbol)
		mut = mut.filtered[unique(cdd.this$mut.idx), ]
		# dim(mut)
		
		mut.agg = aggregate(sample ~ type + mut_type, data=mut, 'length')
		mut.agg = data.frame(name=acc, type=mut.agg$type, mut_type=mut.agg$mut_type, cnt=mut.agg$sample)
		mut.obs = organize.obs(mut.agg)
		expected = results.seq_cat[which(results.seq_cat$acc == acc), 1:4]
		expected[is.na(expected)] = 0
		total.cnt <- colSums(expected[, 2:4]);
		expected <- rbind(expected, data.frame(mut_type='fs_ins', cnt_syn=total.cnt[1], cnt_nonsyn=total.cnt[2], cnt_nonsense=total.cnt[3]));
		expected <- rbind(expected, data.frame(mut_type='fs_del', cnt_syn=total.cnt[1], cnt_nonsyn=total.cnt[2], cnt_nonsense=total.cnt[3]));
		expected <- rbind(expected, data.frame(mut_type='inf_ins', cnt_syn=total.cnt[1], cnt_nonsyn=total.cnt[2], cnt_nonsense=total.cnt[3]));
		expected <- rbind(expected, data.frame(mut_type='inf_del', cnt_syn=total.cnt[1], cnt_nonsyn=total.cnt[2], cnt_nonsense=total.cnt[3]));
		obs_exp <- merge(mut.obs, expected, by='mut_type', all=T)
		obs_exp <- condense.type(obs_exp);
		colnames(obs_exp) <- c('mut_type', 'syn.obs', 'missense.obs', 'nonsense.obs', 'fs.obs', 'inframe.obs', 'syn.exp', 'missense.exp', 'nonsense.exp');
		opt <- compute.selection(obs_exp, joint=T, indel=T)[c(1:6, 8, 10)];	
		result.opt = rbind(result.opt, opt)

		## calculate empirical p-value
		N = N.sim
		opt.perm.check = c(i, rep(NA, 8))
		if(sum(expected[, 2:4]) > 0) {
			if(sum(obs_exp[, c('missense.obs', 'nonsense.obs', 'fs.obs')]) >= 10 & (opt[6] > 1 | opt[7] > 1 | opt[8] > 2)) {
				sizes = rowSums(obs_exp[, grep('obs$', colnames(obs_exp))])
				subtotals = rowSums(obs_exp[, grep('exp$', colnames(obs_exp))])
				perm.list = list()
				for(k in 1:(nrow(obs_exp)-4)) {
					if(subtotals[k] > 0) {
						perm = rmultinom(n=N, size=sizes[k], prob=obs_exp[k, c('syn.exp', 'missense.exp', 'nonsense.exp')]/subtotals[k])
					} else {
						perm = matrix(0, nrow=3, ncol=N)
					}
					perm.list[[k]] = perm
				}
				perm.df = do.call(rbind, perm.list)
				cnt.indel = sum(obs_exp$fs.obs) + sum(obs_exp$inframe.obs)
				perm.indel = rbinom(n=N, size=cnt.indel, prob=c(2/3, 1/3))
				cnt.syn = sum(obs_exp$syn.obs)
				cnt.inf = sum(obs_exp$inframe.obs)
				perm.kept.list = list()
				for(k in 1:N) {
					obs.perm = matrix(perm.df[, k], ncol=3, byrow=T)
					obs.perm = rbind(obs.perm, 0, 0, 0, 0)
					obs.perm = cbind(obs.perm, 0, 0)
					obs.perm[nrow(obs.perm)-3, 4] = perm.indel[k] 
					obs.perm[nrow(obs.perm)-1, 5] = cnt.indel - perm.indel[k] 
					if((sum(obs.perm[, 1]) > 0 & sum(obs.perm[, 1]) <= cnt.syn) | (sum(obs_exp$inframe.obs) > 0 & sum(obs.perm[, 5]) < cnt.inf & sum(obs.perm[, 4]) > 3*sum(obs.perm[, 5]))) 
						perm.kept.list[[length(perm.kept.list)+1]] = obs.perm
				}
				opt.perm.check = c(i, rep(-100, 8))
				if(length(perm.kept.list) > 0) {
					if(length(perm.kept.list) > 100) cat(i, '...', length(perm.kept.list), '\n'); flush.console();
					opt.perm = c()
					for(k in 1:length(perm.kept.list)) {
						obs.perm = perm.kept.list[[k]]
						obs.perm = cbind(obs_exp[, 1], obs.perm, obs_exp[, grep('exp$', colnames(obs_exp))])
						colnames(obs.perm) = colnames(obs_exp)
						rownames(obs.perm) = rownames(obs_exp)
						s = compute.selection(obs.perm, joint=T, indel=T)[c(1:6, 8, 10)];
						opt.perm = rbind(opt.perm, c(i, s))
					}
					idx.mis = which(abs(opt.perm[, 6]) > abs(opt[6]))
					idx.non = which(abs(opt.perm[, 7]) > abs(opt[7]))
					idx.both = which(abs(opt.perm[, 8]) > abs(opt[8]))
					# idx.mis = which(opt.perm[, 6] > opt[6])
					# idx.non = which(opt.perm[, 7] > opt[7])
					# idx.both = which(opt.perm[, 8] > opt[8])
					opt.perm.check = opt.perm[unique(c(idx.mis, idx.non, idx.both)), ]
				}
			}
		}
		result.opt.perm[[i]] = opt.perm.check	
	}
	cat('\n');
	colnames(result.opt) = c('syn.obs', 'missense.obs', 'nonsense.obs', 'fs.obs', 'inframe.obs', 'sel.missense', 'sel.nonsense', 'sel.both')
	# dim(result.opt); length(result.opt.perm)  ## 8342
	cnt.fp = c()
	for(i in 1:length(result.opt.perm)) { 
		perm = result.opt.perm[[i]]
		if(length(perm) == 9) 	perm = matrix(perm, ncol=9, byrow=T)
		fp.mis = fp.non = fp.both = NA
		if(length(which(is.na(perm))) < 5) {
			perm = data.frame(perm)
			colnames(perm) = c('i', colnames(result.opt))
			if(nrow(perm) == 0 | is.na(perm[1, 'sel.missense']) | perm[1, 'sel.missense'] == -100) {
				fp.mis = fp.non = fp.both = 0
			} else {
				fp.mis = length(which(result.opt[i, 6]*(perm$sel.missense-result.opt[i, 6]) > 0))
				fp.non = length(which(result.opt[i, 7]*(perm$sel.nonsense-result.opt[i, 7]) > 0))
				fp.both = length(which(result.opt[i, 8]*(perm$sel.both-result.opt[i, 8]) > 0))
			}
		}
		cnt.fp = rbind(cnt.fp, c(fp.mis, fp.non, fp.both))
	}
	# dim(cnt.fp)  ## 8342    3
	colnames(cnt.fp) = c('fp.missense', 'fp.nonsense', 'fp.both')
	result.opt.fp = data.frame(db_xref=cdd.uniq, cbind(result.opt, cnt.fp), stringsAsFactors=F)
	cat('   ', length(which(!is.na(result.opt.fp$sel.both))), ' regions have selection coefficients estimated.\n');  ## 8342
	
	write.table(result.opt.fp, selection.region.file.name, sep='\t', row.names=F, quote=F)
}

#' Filter the CDD regions for those under significant positive selection. It requires the output.folder contains files from running the parse.aggregated.mut(), annotate.by.cdd.region(), and compute.selection.by.region() functions. That is, prefix.mut.filtered.txt, prefix.mut.region.txt, and prefix.selection.region.txt files must exist.
#' @param output.folder: the folder to write output files
#' @param output.prefix: prefix used to name the output file as prefix.sig.region.txt. 
#' @param N.sim: the number of simulations performed when calling the compute.selection.by.region() function. Default value is 1,000.
#' @param fp.cutoff.rate: false positive cutoff rate (default: 0.1).
#' @param cnt.mis.cutoff: the minimum number of tumors harboring missense mutations in a CDD region (default: 6).
#' @param cnt.non.cutoff: the minimum number of tumors harboring proteinitruncating mutations in a CDD region (default: 6).
#' @param cnt.mis.cutoff: the minimum number of tumors harboring missense or protein-truncating mutations in a CDD region (default: 6).
#' @param cnt.cutoff.rate: the minimum fraction of tumors harboring mutations in the CDD region (default: NULL). If set (range: 0 - 1), the minimum number of tumors will be calculated as round(cnt.cutoff.rate * number of samples in the prefix.mut.filtered.txt file. If this number is larger than cnt.both.cutoff, it will overwrite cnt.both.cutoff, cnt.mis.cutoff, and cnt.non.cutoff. If not set (NULL), cnt.mis.cutoff, cnt.non.cutoff, and cnt.both.cutoff must be set.  
#' @param sel.cutoff: the minimum selection coefficient of a CDD region (default: 1).
#' @return NULL
#' @examples find.sig.selection(output.folder='./examples/', output.prefix='TCGA.CHOL')
#' @export
find.sig.selection = function(output.folder=output.folder, output.prefix=output.prefix, N.sim=1000, fp.cutoff.rate=0.1, cnt.cutoff.rate=NULL, cnt.mis.cutoff=6, cnt.non.cutoff=6, cnt.both.cutoff=6, sel.cutoff=1) {
	cat('filtering CDD regions for those with significant selection coefficients ...\n'); flush.console();
	working.dir = paste(output.folder, '/', sep='');
	mut.file.name = paste(working.dir, output.prefix, '.mut.filtered.txt', sep='');
	selection.region.file.name = paste(working.dir, output.prefix, '.selection.region.txt', sep='');
	sig.region.file.name = paste(working.dir, output.prefix, '.sig.region.txt', sep='');

	mut.filtered = read.table(mut.file.name, sep='\t', header=T, stringsAsFactors=F)
	n.tumors = length(unique(mut.filtered$sample))
	cat('   total tumor count', n.tumors, '\n');
	if(!is.null(cnt.cutoff.rate)) {
		cnt.cutoff = round(n.tumors * cnt.cutoff.rate)
		if(cnt.cutoff > cnt.both.cutoff) cnt.mis.cutoff=cnt.non.cutoff=cnt.both.cutoff=cnt.cutoff
	}
	cat('   cnt.cutoff is', cnt.mis.cutoff, '(missense)', cnt.non.cutoff, '(nonsense)', cnt.both.cutoff, '(both)\n'); flush.console()

	result.opt.fp = read.table(selection.region.file.name, sep='\t', header=T, stringsAsFactors=F);

	fp.cutoff = N.sim * fp.cutoff.rate
	cat('   fp.cntoff is', fp.cutoff, 'out of', N.sim, '\n')
	idx.mis = which(!is.na(result.opt.fp$sel.missense) & abs(result.opt.fp$sel.missense) >= sel.cutoff & result.opt.fp$fp.missense <= fp.cutoff & result.opt.fp$missense.obs >= cnt.mis.cutoff)
	idx.non = which(!is.na(result.opt.fp$sel.nonsense) & abs(result.opt.fp$sel.nonsense) >= sel.cutoff & result.opt.fp$fp.nonsense <= fp.cutoff & (result.opt.fp$nonsense.obs + result.opt.fp$fs.obs) >= cnt.non.cutoff)
	idx.both = which(!is.na(result.opt.fp$sel.both) & abs(result.opt.fp$sel.both) >= sel.cutoff & result.opt.fp$fp.both <= fp.cutoff & (result.opt.fp$missense.obs + result.opt.fp$nonsense.obs + result.opt.fp$fs.obs) >= cnt.both.cutoff)
	result.opt.fp[, 'sig.missense'] = 0; result.opt.fp[idx.mis, 'sig.missense'] = 1;
	result.opt.fp[, 'sig.nonsense'] = 0; result.opt.fp[idx.non, 'sig.nonsense'] = 1;
	result.opt.fp[, 'sig.both'] = 0; result.opt.fp[idx.both, 'sig.both'] = 1;
	idx = unique(c(idx.mis, idx.non, idx.both))
	cat('     ', length(idx), 'passing filters.\n');  ## 35
	sig = cbind(idx=idx, result.opt.fp[idx, ])
	result = sig		
		
	result$pssm = as.numeric(gsub('CDD:', '', result$db_xref))
	result.anno.g = merge(result, unique(anno.g[, c('pssm', 'cdd', 'cdd.name')]), by='pssm', all.x=T, all.y=F)
	result.anno = merge(result.anno.g, unique(anno[, c('pssm', 'cdd', 'cdd.name')]), by='pssm', all.x=T, all.y=F)
	cat('     ', nrow(result.anno), 'after adding CDD anno.g annotations.\n');  ## 81
	result.anno = merge(result.anno, unique(cdd.region[, c('db_xref', 'note')]), by='db_xref', all.x=T, all.y=F)
	cat('     ', nrow(result.anno), 'after adding CDD genbank note annotations.\n');  ## 81
	
	write.table(result.anno, sig.region.file.name, sep='\t', row.names=F, quote=F)
}

#' Plot mutation distribution in CDD regions and save to a PDF file. It requires the output.folder contains files from running the parse.aggregated.mut(), annotate.by.cdd.region(), and find.sig.selection() functions. That is, prefix.mut.filtered.txt, prefix.mut.region.txt, and prefix.sig.region.txt files must exist.
#' @param pdf.file.name: name of the PDF file for saving the plot (default: 'plot.regions.sig.pdf')
#' @param output.folder: the folder to write output files
#' @param output.prefix: prefix used to name the prefix.mut.filtered.txt and prefix.mut.region.txt files. 
#' @param selected.regions: a vector of CDD IDs (db_xref) to plot (default NULL). If set to NULL, all CDD regions in the prefix.sig.region.txt file will be plotted. If it is not NULL, those overalpping with the prefix.sig.region.txt file will be plotted.
#' @param rare.mut.cutoff: the minimum number of missense or protein-truncating mutations (default: 3) in a gene for which to be considered as "rarely" mutated. For a rarely mutated gene, its gene name will be enclosed inside a parenthesis. 
#' @param plot.nrow and plot.col: the number of rows and columns in the plot (default: 3x4 grid)
#' @param main: strings (default: '') appended to the title of the plot, which is the CDD ID.
#' @export
plot_regions = function(output.folder, output.prefix, selected.regions=NULL, pdf.file.name='plot.regions.sig.pdf', rare.mut.cutoff=3, plot.nrow=3, plot.ncol=4, main='') {
	cat('plotting CDD regions ...\n'); flush.console();
	working.dir <- paste(output.folder, '/', sep='');
	mut.file.name <- paste(working.dir, output.prefix, '.mut.filtered.txt', sep='');
	mut.region.file.name <- paste(working.dir, output.prefix, '.mut.region.txt', sep='');
	sig.region.file.name <- paste(working.dir, output.prefix, '.sig.region.txt', sep='');
	if(!is.null(pdf.file.name)) pdf.file.name = paste(working.dir, pdf.file.name, sep='') 

	mut.filtered = read.table(mut.file.name, sep='\t', header=T, stringsAsFactors=F)
	mut.cdd <- read.table(mut.region.file.name, sep='\t', quote='', header=T, stringsAsFactors=F);	
	sig <- read.table(sig.region.file.name, sep='\t', quote='', header=T, stringsAsFactors=F);
	if(!is.null(selected.regions)) {
		sig = sig[which(sig$db_xref %in% selected.regions), ]
	}	
	cat('   ', nrow(sig), 'regions to plot.\n');
	if(nrow(sig) > 0) {
		if(!is.null(pdf.file.name)) pdf(pdf.file.name, width=20, heigh=12)
		if(!is.null(plot.nrow)) par(mfrow=c(plot.nrow, plot.ncol))
		acc = unique(sig$db_xref)  ## some CDDs are repeated because of different values in the "note" column.
		uu = names(seq.gapped.list)
		length(acc) 
		cdd.region.target = c()
		for(target in acc) { 
			one.cdd.region = plot_this_region(target, mut.filtered, mut.cdd, uu, sig, rare.mut.cutoff=rare.mut.cutoff, main=main);
			cdd.region.target = rbind(cdd.region.target, one.cdd.region)
		}
		if(!is.null(pdf.file.name)) dev.off()
		failed = unique(cdd.region.target[which(cdd.region.target$refseq == ''), 'db_xref'])
		cat(length(failed), ' regions failed to plot\n  ', paste(failed, collapse=' '), '\n'); flush.console();
	}
}

#' Plot mutation distribution for one CDD region.
#' @param target: the ID (db_xref) of the CDD region to be plotted.
#' @param mut.filtered: the data frame from the prefix.mut.fitlered.txt file. 
#' @param mut.cdd: the data frame from the prefix.mut.region.txt file. 
#' @param uu: the vector of unique CDD IDs from the seq.gapped.list. 
#' @param sig: the data frame from the prefix.sig.region.txt file
#' @param rare.mut.cutoff: the minimum number of missense or protein-truncating mutations (default: 3) in a gene for which to be considered as "rarely" mutated. For a rarely mutated gene, its gene name will be enclosed inside a parenthesis. 
#' @param main: strings (default: '') appended to the title of the plot, which is the CDD ID.
#' @param verbose: TRUE or FALSE (default: ALSE) controls if the name of the genes in a CDD region will be printed in the standard output.
#' @export
plot_this_region = function(target, mut.filtered, mut.cdd, uu, sig, rare.mut.cutoff=3, main='', verbose=F) {
	plot.drawn = F
	cdd.this = mut.cdd[which(mut.cdd$db_xref == target), ]
	cdd.this = unique(cdd.this[, c('refseq', 'symbol', 'seq.aa', 'mut.aa.coord', 'mut.idx', 'note')])  ## still have redundant rows of the same mut.idx  ## only the first one is take during plotting.
	dim(cdd.this);  length(unique(cdd.this$mut.idx)); 
	cdd.region.target = data.frame(symbol='', refseq='', start=0, end=0, stringsAsFactors=F)
	
	if(nrow(cdd.this) > 0) {
		ss.all = unique(cdd.this$symbol)
		if(nrow(cdd.this) > 2000) cdd.this = cdd.this[sample(1:nrow(cdd.this), 2000), ]
		cdd.region.target = data.frame(symbol=unique(cdd.this$symbol), refseq='', start=0, end=0, stringsAsFactors=F)
		if(length(which(uu == target)) == 1) {  ## some CDD is not included in the sequence alignments.
			target.aligned = seq.gapped.list[[which(uu == target)]] 
			mut.coord.adj = c()
			if(length(target.aligned) > 2) {  ## make sure the target.aligned has sequence alignments
				consensus.aa = target.aligned[length(target.aligned)/2 + 1];
				if(nchar(consensus.aa) == 0) {  ## sometimes the consensus sequence is empty
					seq.idx = which(nchar(target.aligned) > 0)
					seq.idx = seq.idx[which(seq.idx > length(target.aligned)/2)]  ## sometimes the consensus sequence is empty
					if(length(seq.idx) > 0) {
						consensus.aa = target.aligned[min(seq.idx)]
					}
				}
				temp.cnt = 0
				cdd.this.list = split(cdd.this, cdd.this$refseq)
				for(i in 1:length(cdd.this.list)) {
					cdd.this.sub = cdd.this.list[[i]]
					cdd.this.aa = cdd.this.sub[1, 'seq.aa']
					cdd.this.refseq = cdd.this.sub[1, 'refseq']
					refseq.idx = grep(paste0('\\|', cdd.this.refseq, '\\|'), target.aligned)
					refseq.aligned = c()
					if(length(refseq.idx) > 0) {
						refseq.aligned = target.aligned[length(target.aligned)/2 + refseq.idx]
					} else {  ## I decided not to do the alignment because those unaligned refseqs may have an aligned refseq representing the same gene. 
						temp.cnt = temp.cnt + 1
						if(temp.cnt < 50) {
							# cat(temp.cnt, '...\n'); flush.console();
							align = pwalign::pairwiseAlignment(pattern=cdd.this.aa, subject=consensus.aa, type='global')
							consensus.aligned = pwalign::toString(pwalign::alignedSubject(align))
							if(consensus.aligned == consensus.aa) {
								refseq.aligned = pwalign::toString(pwalign::alignedPattern(align))
								# cat(target, 'vs.', cdd.this.refseq, 'alignment succeeded\n'); flush.console();
							}
						}
					}				
					if(length(refseq.aligned) > 0) {			
						for(j in 1:nrow(cdd.this.sub)) {
							one = cdd.this.sub[j, ]
							gaps = gregexpr('-', refseq.aligned)[[1]]
							one.coord = one$mut.aa.coord + length(gaps[which(gaps < one$mut.aa.coord)]) 
							one$mut.coord.adj = one.coord
							mut.coord.adj = rbind(mut.coord.adj, one)
						}
					} else {
						cdd.this.sub$mut.coord.adj = NA
						mut.coord.adj = rbind(mut.coord.adj, cdd.this.sub)				
					}
				}
			} else {
				mut.coord.adj = cdd.this
				mut.coord.adj$mut.coord.adj = NA
			}
			
			if(length(which(!is.na(mut.coord.adj$mut.coord.adj))) == 0) {
				agg = aggregate(mut.idx ~ refseq, data=cdd.this, 'length')
				represent = agg[which(agg$mut.idx == max(agg$mut.idx)), 'refseq']
				np.idx = grep('NP_', represent)
				if(length(np.idx) == 0) {
					represent = represent[which.min(nchar(represent))]
				} else if(length(np.idx) == 1) {
					represent = represent[np.idx]
				} else if(length(np.idx) > 1) {
					represent = represent[np.idx]
					represent = represent[which.min(nchar(represent))]
				}
				mut.coord.adj = cdd.this[which(cdd.this$refseq == represent), ]
				mut.coord.adj$mut.coord.adj = mut.coord.adj$mut.aa.coord
				cat('Warning: no alignment for', target, '. A representative protein is used.\n'); flush.console()
				main = paste('#', length(ss.all), '#', main)
			}
				
			if(length(which(!is.na(mut.coord.adj$mut.coord.adj))) > 0) {
				cdd.region.target = unique(cdd.region[which(cdd.region$db_xref == target & cdd.region$refseq %in% cdd.this$refseq), c('symbol', 'refseq', 'start', 'end')])
				cdd.this = merge(cdd.this, mut.coord.adj[, c('refseq', 'mut.idx', 'mut.coord.adj')], by=c('refseq', 'mut.idx'), all=F)
				mut.adj = data.frame(coord.adj=cdd.this$mut.coord.adj, mut.idx=cdd.this$mut.idx, type=mut.filtered[cdd.this$mut.idx, 'type'], symbol=mut.filtered[cdd.this$mut.idx, 'symbol'], refseq=cdd.this$refseq, stringsAsFactors=F)
				mut.adj$color = ifelse(mut.adj$type == 'syn', 'green', ifelse(mut.adj$type %in% c('inf_ins', 'inf_del'), 'cyan3', ifelse(mut.adj$type == 'missense', 'blue', 'red')))
				## remove redundant rows for the same mut.idx
				a = split(mut.adj, mut.adj$mut.idx)
				mut.adj = do.call(rbind, lapply(a, function(x) x[1, ]))
				tb = as.data.frame(table(mut.adj$coord.adj, mut.adj$color), stringsAsFactors=F)
				tb$Var1 = as.numeric(tb$Var1)
				tb$Var2 = factor(tb$Var2, levels=c('green', 'cyan3', 'blue', 'red'))
				tb = tb[order(tb$Var2), ]
				tb$Var2 = as.character(tb$Var2)
				tb$pch = ifelse(tb$Var2 == 'blue', 19, ifelse(tb$Var2 == 'cyan3', 2, ifelse(tb$Var2 == 'red', 4, 3)))
				cdd.length = max(nchar(target.aligned[(length(target.aligned)/2 + 1):length(target.aligned)]))
				cdd.length = max(cdd.length, max(cdd.region.target$end - cdd.region.target$start))
				plot(tb$Var1, tb$Freq, type='h', col=tb$Var2, xlim=c(1, max(cdd.length, max(tb$Var1) + 5)), ylim=c(0, max(tb$Freq) + 5), main=paste(main, target), xlab='', ylab='Frequency', las=1)
				plot.drawn = T				
				note = cdd.this[1, 'note']; if(nchar(note) > 50) note = paste0(substr(note, 1, 50), '...')
				tb = tb[which(tb$Freq > 0), ]
				points(tb$Var1, tb$Freq, pch=tb$pch, col=tb$Var2)
				mtext(note)
				mtext('Position', side=1, line=2)
				ss.mut = unique(mut.adj$symbol)
				ss.missing = ss.all[which(!ss.all %in% ss.mut)]
				ss = sapply(ss.mut, function(x) if(length(which(mut.adj$symbol == x & mut.adj$type != 'syn')) < rare.mut.cutoff) paste0('(', x,')') else x)
				ss = paste(ss, collapse=' ')
				if(length(ss.missing) > 0) ss = paste0(ss, ' /', paste(ss.missing, collapse=' '), '/') 
				if(length(ss.all) > 5) {
					mtext(paste(c(paste0('[', length(ss.all), ']'), strsplit(ss, ' ')[[1]][1:5], '...'), collapse=' '), side=1, line=3)
				} else {
					mtext(paste(c(paste0('[', length(ss.all), ']'), ss), collapse=' '), side=1, line=3)
				}
				this.result = sig[which(sig$db_xref == target), ]
				if(nrow(this.result) > 0) {
					label.mis = paste0('mis [', round(this.result$sel.missense, 2), ']; '); if(this.result[1, 'sig.missense'] == 1) label.mis = paste('*', label.mis);
					label.non = paste0('non [', round(this.result$sel.nonsense, 2), ']; '); if(this.result[1, 'sig.nonsense'] == 1) label.non = paste('*', label.non);
					label.mis.non = paste0('both [', round(this.result$sel.both, 2), ']'); if(this.result[1, 'sig.both'] == 1) label.mis.non = paste('*', label.mis.non);
					mtext(paste(label.mis, label.non, label.mis.non, sep=''), side=3, line=-2)
				}
				if(verbose) cat(target, '...', ss, '\n')
			
				agg = aggregate(color ~ symbol + refseq, data=mut.adj, 'length')
				selected.refseq = unlist(lapply(split(agg, agg$symbol), function(x) x[which.max(x$color), 'refseq']))
				temp.1 = cdd.region.target[which(cdd.region.target$refseq %in% selected.refseq), ]
				# cdd.region.target = unique(cdd.region[which(cdd.region$db_xref == target & cdd.region$refseq %in% mut.adj$refseq), c('symbol', 'refseq', 'start', 'end')])
				# temp = merge(cdd.region.target, unique(mut.adj[, c('symbol', 'refseq')]), by='refseq', all=T)  ## in case some refseq are not in cdd.region
				# temp = temp[which(temp$refseq %in% selected.refseq), ]
				# cdd.region.target = data.frame(db_xref=target, refseq=temp$refseq, symbol=temp$symbol.y, start=temp$start, end=temp$end, stringsAsFactors=F)
				missing.ss = ss.all[which(!ss.all %in% temp.1$symbol)]
				if(length(missing.ss) > 0) {
					cdd.this = unique(mut.cdd[which(mut.cdd$db_xref == target), c('refseq', 'symbol', 'seq.aa', 'mut.aa.coord', 'mut.idx', 'note')])
					agg = aggregate(mut.idx ~ symbol + refseq, data=cdd.this[which(cdd.this$symbol %in% missing.ss), ], 'length')
					selected.refseq = unlist(lapply(split(agg, agg$symbol), function(x) x[which.max(x$mut.idx), 'refseq']))
					temp.2 = unique(cdd.region[which(cdd.region$refseq %in% selected.refseq & cdd.region$db_xref == target), c('symbol', 'refseq', 'start', 'end')])
					temp.2$refseq = ''
					temp.1 = rbind(temp.1, temp.2)
				}
				cdd.region.target = temp.1			
			}
		}
		
		failed.idx = which(cdd.region.target$refseq == '')
		if(length(failed.idx) == nrow(cdd.region.target)) {
			cat('failed to plot region', target, paste(unique(cdd.this$symbol), collapse=' '), '\n');
		} else if(length(failed.idx) > 0) {
			cat('missing some genes for', target, paste(unique(cdd.region.target[failed.idx, 'symbol']), collapse=' '), '\n');
		}

		cdd.region.target$db_xref = target
		
	}
	if(!plot.drawn) {
		plot(0, type='n', xlim=c(0, 1), ylim=c(0, 1), xlab='', ylab='', xaxt='n', yaxt='n', main=paste(main, target))
		text(0.5, 0.5, 'No mutation distribution is available for this domain.')
	}
	return(cdd.region.target)
}

#' Perform DUST analysis to estimate selection coefficients on protein domains and identify those under significant positive selection.
#' @param input.file.name: A VCF-formatted file to read SnpEff annotated somatic variants. This input file shall contain these fields: Tumor_Sample_Barcode, Chromosome, Start_Position, dbSNP_RS, Reference_Allele, Tumor_Seq_Allele2, FILTER, One_Consequence, Hugo_Symbol, Gene, Feature, ENSP, HGVSc, HGVSp_Short, Amino_acids, Codons, ENSP, RefSeq, Entrez_Gene_Id.
#' @param output.folder: the folder to write the temporary and final result files. 
#' @param output.prefix: prefix used to name the temporary and final prediction files. Temporary and output files include prefix.outlier.txt, prefix.mut.filtered.txt, prefix.mut.region.txt, prefix.error.txt, prefix.exp_seq.by_type.region.txt, prefix.selection.region.txt, prefix.sig.region.txt, and a plot pdf file.
#' @param steps: A vector of integers indicating which functions to execute. 1-find.outliers(), 2-parse.aggregated.mut(), 3-annotate.by.cdd.region(), 4-calculate.expected.by.region, 5-compute.selection.by.region(), 6-find.sig.selection(), 7-plot.regions(). Default to 1:7.
#' @param N.sim: the number of simulations performed when calling the compute.selection.by.region() function. Default value is 1,000.
#' @param fp.cutoff.rate: false positive cutoff rate (default: 0.1).
#' @param cnt.mis.cutoff: the minimum number of tumors harboring missense mutations in a CDD region (default: 6).
#' @param cnt.non.cutoff: the minimum number of tumors harboring proteinitruncating mutations in a CDD region (default: 6).
#' @param cnt.mis.cutoff: the minimum number of tumors harboring missense or protein-truncating mutations in a CDD region (default: 6).
#' @param cnt.cutoff.rate: the minimum fraction of tumors harboring mutations in the CDD region (default: NULL). If set (range: 0 - 1), the minimum number of tumors will be calculated as round(cnt.cutoff.rate * number of samples in the prefix.mut.filtered.txt file. If this number is larger than cnt.both.cutoff, it will overwrite cnt.both.cutoff, cnt.mis.cutoff, and cnt.non.cutoff. If not set (NULL), cnt.mis.cutoff, cnt.non.cutoff, and cnt.both.cutoff must be set.  
#' @param rare.mut.cutoff: the minimum number of missense or protein-truncating mutations (default: 3) in a gene for which to be considered as "rarely" mutated. For a rarely mutated gene, its gene name will be enclosed inside a parenthesis. 
#' @return NULL
#' @examples gust(input.file.name='./examples/TCGA.CHOL.somatic.maf.gz', output.folder='./examples/', output.prefix='TCGA.CHOL');
#' @export
dust <- function(input.file.name, output.folder, output.prefix, steps=1:7, N.sim=1000, fp.cutoff.rate=0.1, cnt.cutoff.rate=NULL, cnt.mis.cutoff=6, cnt.non.cutoff=6, cnt.both.cutoff=6, rare.mut.cutoff=3) {
	output.folder <- paste(output.folder, '/', sep='');
	if(1 %in% steps) {
		find.outliers(input.file.name=input.file.name, output.folder=output.folder, output.prefix=output.prefix);
	}
	if(2 %in% steps) {
		parse.aggregated.mut.2024(input.file.name=input.file.name, output.folder=output.folder, output.prefix=output.prefix, fa_aa=fasta_aa, fa_bp=fasta_cdf, recover=F);
	}
	if(3 %in% steps) {
		annotate.by.cdd.region(output.folder=output.folder, output.prefix=output.prefix) 
	}
	if(4 %in% steps) {
		calculate.expected.by.region(output.folder=output.folder, output.prefix=output.prefix)
	}
	if(5 %in% steps) {
		compute.selection.by.region(output.folder=output.folder, output.prefix=output.prefix, N.sim=N.sim)
	}
	if(6 %in% steps) {
		find.sig.selection(output.folder=output.folder, output.prefix=output.prefix, N.sim=N.sim, fp.cutoff.rate=fp.cutoff.rate, cnt.cutoff.rate=cnt.cutoff.rate, cnt.mis.cutoff=cnt.mis.cutoff, cnt.non.cutoff=cnt.non.cutoff, cnt.both.cutoff=cnt.both.cutoff)
	}
	if(7 %in% steps) {
		fn = paste0(output.prefix, '.sig.region.txt')
		selected = c()
		if(file.exists(paste0(output.folder, '/', fn))) selected = read.table(paste0(output.folder, '/', fn), sep='\t', quote='', header=T, stringsAsFactors=F)
		if(!is.null(selected) && nrow(selected) > 0) {
			plot_regions(output.folder=output.folder, output.prefix=output.prefix, selected.regions=selected$db_xref, pdf.file.name=gsub('.txt', 'plot.pdf', fn)) ## some CDD has fewer muts plotted than muts observed because the refseq is not found in the alignments.
		}
	}
}



