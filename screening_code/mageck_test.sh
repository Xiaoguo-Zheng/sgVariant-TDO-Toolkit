#! /bin/sh
awk -F "\t" '{$3="";print $0}' Final_sgRNA_Count_Matrix.xls |sed 's/\t\t/\t/g' >Table1.txt

mageck test -k Table1.txt \
			-n H1H2_vs_W1W2  \
			-t H1,H2 -c W1,W2 \
			--norm-method total \
			--adjust-method fdr \
			--sort-criteria neg \
			--remove-zero both \
			--remove-zero-threshold 0 \
			--pdf-report \
			--normcounts-to-file 	\
			--keep-tmp
			
mageck test -k Table1.txt \
			-n M1M2_vs_W1W2  \
			-t M1,M2 \
			-c W1,W2 \
			--norm-method total \
			--adjust-method fdr \
			--sort-criteria neg \
			--remove-zero both \
			--remove-zero-threshold 0 \
			--pdf-report \
			--normcounts-to-file 	\
			--keep-tmp

mageck test -k Table1.txt \
			-n L1L2_vs_W1W2  \
			-t L1,L2 \
			-c W1,W2 \
			--norm-method total \
			--adjust-method fdr \
			--sort-criteria neg \
			--remove-zero both \
			--remove-zero-threshold 0 \
			--pdf-report \
			--normcounts-to-file 	\
			--keep-tmp
			
mageck test -k Table1.txt \
			-n N1N2_vs_W1W2  \
			-t N1,N2 \
			-c W1,W2 \
			--norm-method total \
			--adjust-method fdr \
			--sort-criteria neg \
			--remove-zero both \
			--remove-zero-threshold 0 \
			--pdf-report \
			--normcounts-to-file 	\
			--keep-tmp			
			
						
mageck test -k Table1.txt \
			-n H1H2_vs_N1N2  \
			-t H1,H2 \
			-c N1,N2 \
			--norm-method total \
			--adjust-method fdr \
			--sort-criteria neg \
			--remove-zero both \
			--remove-zero-threshold 0 \
			--pdf-report \
			--normcounts-to-file 	\
			--keep-tmp			
			