!Sample_data_processing = Adapters were trimmed off from raw reads with Trimmomatic with argument "ILLUMINACLIP:$FA_ADAPTER:6:30:10 LEADING:3 TRAILING:3 MINLEN:36 SLIDINGWINDOW:4:15". 
!Sample_data_processing = Raw reads were aligned with Hisat2 with arguments "--no-mixed --rna-strandness RF --dta --fr" to produce a SAM file.
!Sample_data_processing = Duplicate reads were removed with Picard using default
setting
!Sample_data_processing = Alignments in SAM file were assembled into transcripts abundances with stringtie with argument "--rf".


!Sample_processed_data_files_format_and_content = .stringtie.count: TSV table containing abundance of transcripts. 
!Sample_supplementary_file = {processedFile_indexed}

