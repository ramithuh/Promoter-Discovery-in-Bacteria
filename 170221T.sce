// Clear Variables and Load helper functions
clear;
exec('helpers.sce'); 

function [extended_strand] = extend_strand(strand)
    // For each gene in the protein table array modify its start and end indexes such that
    // we obtain 50 bases upstream and 3 bases downstream

    for i=1:length(strand)/2
        strand(i,1) = strand(i,1) - 50;
        strand(i,2) = strand(i,2) + 3;
    end
    extended_strand = strand;

endfunction

function [ filtered_genes, ignore_rows ] = ignore_genes(strand, verbose )
    if ~exists("verbose","local") then
        verbose = 0
    end
    

    //A gene is ignored if the previous gene is less than 50 bases upstream of it
    ignored_count = 0
    overlap_count = 0

    ignore_rows = [];
    filtered_genes = [] //variable to store the filtered genes
    
    i = 1;

    while(i <= length(strand)/2) //loop through each gene
        if(i==1)
            filtered_genes = [filtered_genes;strand(i,1),strand(i,2)];
        else
            last_pos      = length(filtered_genes)/2;
            last_gene_end = filtered_genes(last_pos,2) //'end' of the last_gene that was added to the filtered genes

            if(strand(i,1) - last_gene_end -1 < 50) //closer than 50 bases
                ignored_count = ignored_count+1;
                ignore_rows = [ignore_rows; strand(i,1),strand(i,2)]; //save ignored start,end positions
                    
                if(strand(i,1) - last_gene_end < 0)  //overlapping genes
                    overlap_count = overlap_count+1;
                end

            else //if the gap between the genes is greater than 50 bases
                filtered_genes = [filtered_genes;strand(i,1),strand(i,2)];
            end
        end

        i = i+1;
    end

    if(verbose)
        disp('Total Genes : ' + string(length(strand)/2));
        disp('Ignored genes : ' + string(ignored_count));
        disp('  Overlapping genes : ' + string(overlap_count));
        disp('Remaining genes after filtering : ' + string(length(filtered_genes)/2));
    end

endfunction


function [distribution, proportion_ignored] = Question_1(filtered_strand, ignored_genes, fasta_in, strand_type, verbose)
    // Question 1 : 
    // Obtain the distribution of bases between genes and the proportion 
    // of genes that can be neglected for being less than 50 bases downstream 
    // of the next one.

    if(~exists("verbose","local")) then
        verbose = 0;
    end

    total_genes = length(filtered_strand)/2 + length(ignored_genes)/2;

    proportion_ignored  = (length(ignored_genes)/2) /total_genes;

    bases = []

    //loop through the filtered genes and update the bases array
    for i=1:length(filtered_strand)/2
        b = get_fasta_at(fasta_in, filtered_strand(i,1), filtered_strand(i,2), strand_type);
        bases = [bases, b];
    end

    distribution = get_base_hist(bases);
    if(verbose)
        if(strand_type == 1) then
            strand_type = 'sense'
        else
            strand_type = 'antisense'
        end

        disp('Question 1')
        disp(' ⟶  Distribution of bases in the filtered genes' + ' of the ' + strand_type + ' strand');
        disp(' ⟶ ' + string(distribution));
        disp(' ⟶  Total Genes : ' + string(total_genes));
        disp(' ⟶  Ignored genes : ' + string(length(ignored_genes)/2));
        disp(' ⟶  Proportion of ignored genes ' + string(proportion_ignored));
    end 
endfunction

function [PPM] = Question_2(filtered_strand, fasta_in, verbose)
    // Question 2
    // Perform a W matching local search (for an intact query) to 
    // locate a Pribnow box promoter within upstream positions 5 to 30 
    // of each sequence. 
    // Using the first 1000 sequences, obtain a position probability matrix (PPM) 
    // with 10 positions for the Pribnow box.
    fs = filtered_strand
    u_row = size(filtered_strand,1);

    // Get region with Pribnow Box
    pribnow_start = 30;  // Start position offset of pribnow Box (upstream positions upto 30)
    pribnow_stop  = 5;   // Stop position offset of pribnow Box (upstream positions upto 5)

    pribnow_len = pribnow_start - pribnow_stop; // Length of region with pribnow Box
    pribnow_mat_len = 10; // Length of PPM of pribnow Box

    k_init = 0.01; // Initial probability (to allow log calculations)

    pribnow_pos_freq_matrix = k_init*ones(4,pribnow_mat_len); // Create empty frequency matrix for nucleotides to subsequently generate PPM of pribnow Box

    pribnow_query = ascii('TATAAT'); // pribnow box of Escherichia coli 
        
    // Start search    
    cs = 0

    for n_key=1:1000
        // m_check = get_fasta_at(fasta_in,gp(n_key,1),gp(n_key,1)+3,1);
            
        // Get position of pribnow Box
        pribnow_seq = get_fasta_at(fasta_in, fs(n_key,1)-pribnow_start,fs(n_key,1)-pribnow_stop,1); // Potential region to search

        [ax,ay,pribnow_pos] = traceback_prom(pribnow_seq,pribnow_query,1,-1,gap_penalty); // Promoter alignment match A or T (W) with A or T (W)

        if ((pribnow_pos>0)&(pribnow_pos<(pribnow_len-pribnow_mat_len))) then
            // Update the pribnow position frequency matrix
            cs = cs +1
            pribnow_post_prom_seq = pribnow_seq((pribnow_pos+1):length(pribnow_seq));
            pribnow_pos_freq_matrix = update_pos_freq_mat(pribnow_post_prom_seq,pribnow_pos_freq_matrix,pribnow_mat_len);
        end

    end
    disp('Considered ' + string(cs));


    pribnow_ppm = get_ppm(pribnow_pos_freq_matrix); // Generate pribnow Box PPM
    PPM = pribnow_ppm;
    disp('Pribnow Box PPM');
    disp(pribnow_ppm); 

endfunction

// Load helper functions

exec('helpers.sce'); 

// 1. Preprocessing
    // extract the gene positions from the protein table

    // [gp, gn, ncp, ncn] = get_protein_pos_array('proteins_152_747609.txt');
    // fasta_in = 'NZ_CP046291.1.fasta';

    [gp, gn, ncp, ncn] = get_protein_pos_array('other/NZ_CP046280_1_protein_table.txt');
    fasta_in = 'other/NZ_CP046280_1_genome.fasta';
   
    // Filtering and extending of the sense strand
    [f_gp, i_gp] = ignore_genes(gp);
    e_gp         = extend_strand(f_gp);

    // Filtering and extending of the antisense strand
    [f_gn, i_gn] = ignore_genes(gn);
    e_gn         = extend_strand(f_gn);

// Answer of Question 1 : 

    [dist1, _ ] = Question_1(f_gp, i_gp, fasta_in,1,1);
    [dist2, _ ] = Question_1(f_gn, i_gn, fasta_in,0,1);

    distribution = dist1 + dist2;
    total_ignored_proportion = (length(i_gp)/2 + length(i_gn)/2) / (length(gn)/2 + length(gp)/2);

    disp(" Distribution " + string(distribution));
    disp(" Total igored proportion " + string(total_ignored_proportion));


// Answer of Question 2 : 
    PPM = Question_2(f_gp, fasta_in, 1);

    [w,su]=ppm_info(PPM,[0.25 0.25 0.25 0.25]);

    disp("entropy of PPM");
    col_entropy = sum(su,1); // Get the entropy for each column
    disp(col_entropy);
    // Question_2(e_gn, fasta_in, 1);

// Answer of Question 3 : 


//Answer of Question 4 : Perform statistical alignment of the remaining sequences




//display first row of variables gp and e_gp
// disp(gp(1:5,:));

// disp("ignored..")
// disp(length(i_gp)/2);

// disp(e_gp(1:5,:));


// c = []
// c = get_fasta_at(fasta_in,0,70,1)
// disp(c);
// c = [c; c]
// histogram = get_base_hist(c);
// disp(histogram);