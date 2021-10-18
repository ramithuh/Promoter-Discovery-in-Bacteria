
function g=gap_penalty(k)
    // Gap penalty function
    g_alpha = 0;
    g_beta = 2;
    g = -(g_alpha+g_beta*k);
endfunction

function [pos_scores]=stat_align(seq,ppm)
    // Do a statistical alignment
    [m,n]=size(ppm);
    score_vec = [];
    seq_len = length(seq);
    
    // Do statistical alignment
    for k=1:(seq_len-n+1)
        temp_sum = 0;
        for i=1:n
            pos_base_key = get_base_key(seq(k+i-1));
            temp_sum = temp_sum+log(ppm(pos_base_key,i));
        end
        score_vec = [score_vec,temp_sum];
    end
    pos_scores = score_vec;
endfunction

function [pos_scores,col_entropy,col_keys]=stat_align_entropy(seq,ppm,entropy_thresh)
    // Do a statistical alignment using positions with high entropy
    [m,n]=size(ppm);
    score_vec = [];
    seq_len = length(seq);
    
    // Find entropy of PPM for equiprobable baseline and get columns above the threshold
    [w,su]=ppm_info(ppm,[0.25 0.25 0.25 0.25]);
    col_entropy = sum(su,1); // Get the entropy for each column
    col_pos = 1:n; // Column numbers
    col_keys = col_pos(col_entropy>entropy_thresh); // Get keys of columns above the required entropy threshold
    
    col_key_count = length(col_keys);
    col_key_max = max(col_keys); // Ignore all columns of the PPM after the last column with significant entropy
    
    // Find entropy of consensus to normalize
    con_score = 0;
    for i=1:col_key_count
        con_score = con_score+log(max(ppm(:,col_keys(i)))); // Sum maximum entropy of each significant column
    end
    
    // Do statistical alignment
    for k=1:(seq_len-col_key_max+1)
        temp_sum = 0;
        for i=1:col_key_count
            pos_base_key = get_base_key(seq(k+i-1));
            temp_sum = temp_sum+log(ppm(pos_base_key,col_keys(i)));
        end
        score_vec = [score_vec,temp_sum];
    end
    pos_scores = score_vec;
endfunction

function [consensus_score] = get_consensus_score_(PPM, entropy_thresh)
    [m,n]=size(PPM);
    score_vec = [];
    
    // Find entropy of PPM for equiprobable baseline and get columns above the threshold
    [w,su] = ppm_info(PPM,[0.25 0.25 0.25 0.25]);
    col_entropy = sum(su,1); // Get the entropy for each column
    col_pos = 1:n; // Column numbers
    col_keys = col_pos(col_entropy > entropy_thresh); // Get keys of columns above the required entropy threshold
    
    col_key_count = length(col_keys);
    col_key_max = max(col_keys); // Ignore all columns of the PPM after the last column with significant entropy
    
    // Find entropy of consensus to normalize
    consensus_score = 0;
    for i=1:col_key_count
        consensus_score = consensus_score+log(max(PPM(:,col_keys(i)))); // Sum maximum entropy of each significant column
    end

endfunction

function result=prom_base(base_x,base_y,score_m,score_mm)
   // Match mismatch score
   // score_m - match score
   // score_mm - mismatch score
   
   if (base_x==ascii('A'))|(base_x==ascii('T')) then
       base_x = ascii('W');
   end
   
   if (base_y==ascii('A'))|(base_y==ascii('T')) then
       base_y = ascii('W');
   end
   
   if (base_x==base_y) then
       result = score_m;
   else
       result = score_mm;
   end
endfunction

function result=comp_base(base_x,base_y,score_m,score_mm)
   // Match mismatch score
   // score_m - match score
   // score_mm - mismatch score
   if (base_x==base_y) then
       result = score_m;
   else
       result = score_mm;
   end
endfunction

function matrix_result = scorematrix_global(seq_x,seq_y,score_m,score_mm,gap_penalty)
    // Function to perform global search
    len_x = length(seq_x);
    len_y = length(seq_y);
    basic_mat = zeros(len_x+1,len_y+1);
    
    //Initialize with gap penalties
    for k_x=1:(len_x+1)
        basic_mat(k_x,1)=gap_penalty(k_x-1);
    end
    
    for k_y=1:(len_y+1)
        basic_mat(1,k_y)=gap_penalty(k_y-1);
    end
    
    //Recurrance relation
    for k_x=2:(len_x+1)
        for k_y=2:(len_y+1)
            score_match = basic_mat(k_x-1,k_y-1)+comp_base(seq_x(k_x-1),seq_y(k_y-1),score_m,score_mm);
            score_gap_x = basic_mat(k_x,k_y-1)+gap_penalty(1);
            score_gap_y = basic_mat(k_x-1,k_y)+gap_penalty(1);
            basic_mat(k_x,k_y)=max([score_match score_gap_x score_gap_y]);
        end
    end
    matrix_result = basic_mat;
endfunction

function [align_x,align_y] = traceback_global(seq_x,seq_y,score_m,score_mm,gap_penalty)
    // Traceback for global alignment
    score_matrix = scorematrix_global(seq_x,seq_y,score_m,score_mm,gap_penalty);
    
    len_x = length(seq_x);
    len_y = length(seq_y);
    
    base_x = seq_x;
    base_y = seq_y;
    
    out_x = [];
    out_y = [];
    
    k_x = len_x;
    k_y = len_y;
    
    
    // Perform the traceback
    while (k_x>0&k_y>0) then
        k_gx = score_matrix(k_x,k_y+1);
        k_gy = score_matrix(k_x+1,k_y);
        k_diag = score_matrix(k_x,k_y);
        k_c = score_matrix(k_x+1,k_y+1);
        
        if k_c==k_diag+comp_base(seq_x(k_x),seq_y(k_y),score_m,score_mm) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,base_y(k_y)];
            k_x=k_x-1;
            k_y=k_y-1;
        elseif k_c==k_gx+gap_penalty(1) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,45];
            k_x=k_x-1;
        elseif k_c==k_gy+gap_penalty(1) then
            out_x = [out_x,45];
            out_y = [out_y,base_y(k_y)];
            k_y=k_y-1;
        end
    end
    
    // Write the output
    if (k_x>0&k_y==0) then
        pre_x = ascii(base_x(1:k_x));
        pre_y = ascii(45*ones(1,k_x));
        align_x = strcat([pre_x,strrev(ascii(out_x))]);
        align_y = strcat([pre_y,strrev(ascii(out_y))]);
    elseif (k_x==0&k_y>0) then
        pre_y = ascii(base_y(1:k_y));
        pre_x = ascii(45*ones(1,k_y));
        align_x = strcat([pre_x,strrev(ascii(out_x))]);
        align_y = strcat([pre_y,strrev(ascii(out_y))]);       
    else
        align_x = strrev(ascii(out_x));
        align_y = strrev(ascii(out_y));          
    end
endfunction

function matrix_result = scorematrix_local(seq_x,seq_y,score_m,score_mm,gap_penalty)
    //Function to perform local search
    len_x = length(seq_x);
    len_y = length(seq_y);
    basic_mat = zeros(len_x+1,len_y+1);
    
    //No need to initialize because gap penalty will always be negative
    //Recurrance relation
    for k_x=2:(len_x+1)
        for k_y=2:(len_y+1)
            score_match = basic_mat(k_x-1,k_y-1)+comp_base(seq_x(k_x-1),seq_y(k_y-1),score_m,score_mm);
            score_gap_x = basic_mat(k_x,k_y-1)+gap_penalty(1);
            score_gap_y = basic_mat(k_x-1,k_y)+gap_penalty(1);
            basic_mat(k_x,k_y)=max([score_match score_gap_x score_gap_y 0]); // Add the zero to global search
        end
    end
    matrix_result = basic_mat;
endfunction

function [align_x,align_y,start_pos] = traceback_local(seq_x,seq_y,score_m,score_mm,gap_penalty)
    // Traceback for local alignment
    score_matrix = scorematrix_local(seq_x,seq_y,score_m,score_mm,gap_penalty);
    //disp(score_matrix');
    
    pre_x = 0;
    pre_y = 0;
    
    len_x = length(seq_x);
    len_y = length(seq_y);
    
    base_x = seq_x;
    base_y = seq_y;
    
    out_x = [];
    out_y = [];
    
    [m,n]=max(score_matrix);
    k_x = n(1,1)-1;
    k_y = n(1,2)-1;
    
    post_x = ascii(seq_x((k_x+1):len_x));
    post_y = ascii(seq_y((k_y+1):len_y));
    
    // Perform the traceback
    while (k_x>0&k_y>0) then
        k_gx = score_matrix(k_x,k_y+1);
        k_gy = score_matrix(k_x+1,k_y);
        k_diag = score_matrix(k_x,k_y);
        k_c = score_matrix(k_x+1,k_y+1);
        
        if k_c==k_diag+comp_base(seq_x(k_x),seq_y(k_y),score_m,score_mm) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,base_y(k_y)];
            k_x=k_x-1;
            k_y=k_y-1;
        elseif k_c==k_gx+gap_penalty(1) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,45];
            k_x=k_x-1;
        elseif k_c==k_gy+gap_penalty(1) then
            out_x = [out_x,45];
            out_y = [out_y,base_y(k_y)];
            k_y=k_y-1;
        elseif k_c==0;
            break;
        end
    end
    
    // Write the output
    if (k_x>0&k_y==0) then
        pre_x = ascii(base_x(1:k_x));
        pre_y = ascii(45*ones(1,k_x));
        align_x = strcat([pre_x,strrev(ascii(out_x)),post_x]);
        align_y = strcat([pre_y,strrev(ascii(out_y)),post_y]);
    elseif (k_x==0&k_y>0) then
        pre_y = ascii(base_y(1:k_y));
        pre_x = ascii(45*ones(1,k_y));
        align_x = strcat([pre_x,strrev(ascii(out_x)),post_x]);
        align_y = strcat([pre_y,strrev(ascii(out_y)),post_y]);       
    else
        align_x = strrev(ascii(out_x));
        align_y = strrev(ascii(out_y));          
    end
    start_pos = length(pre_y); 
endfunction

function matrix_result = scorematrix_prom(seq_x,seq_y,score_m,score_mm,gap_penalty)
    //Function to perform local search
    len_x = length(seq_x);
    len_y = length(seq_y);
    basic_mat = zeros(len_x+1,len_y+1);
    
    //No need to initialize because gap penalty will always be negative
    //Recurrance relation
    for k_x=2:(len_x+1)
        for k_y=2:(len_y+1)
            score_match = basic_mat(k_x-1,k_y-1)+prom_base(seq_x(k_x-1),seq_y(k_y-1),score_m,score_mm);
            score_gap_x = basic_mat(k_x,k_y-1)+gap_penalty(1);
            score_gap_y = basic_mat(k_x-1,k_y)+gap_penalty(1);
            basic_mat(k_x,k_y)=max([score_match score_gap_x score_gap_y 0]); // Add the zero to global search
        end
    end
    matrix_result = basic_mat;
endfunction

function [align_x,align_y,prom_pos] = traceback_prom(seq_x,seq_y,score_m,score_mm,gap_penalty)
    // Traceback for local alignment
    score_matrix = scorematrix_prom(seq_x,seq_y,score_m,score_mm,gap_penalty);
    pre_x = 0;
    pre_y = 0;
    
    //disp(score_matrix');
    
    len_x = length(seq_x);
    len_y = length(seq_y);
    
    base_x = seq_x;
    base_y = seq_y;
    
    out_x = [];
    out_y = [];
    
    [m,n]=max(score_matrix);
    k_x = n(1,1)-1;
    k_y = n(1,2)-1;
    
    post_x = ascii(seq_x((k_x+1):len_x));
    post_y = ascii(seq_y((k_y+1):len_y));
    
    // Perform the traceback
    while (k_x>0&k_y>0) then
        k_gx = score_matrix(k_x,k_y+1);
        k_gy = score_matrix(k_x+1,k_y);
        k_diag = score_matrix(k_x,k_y);
        k_c = score_matrix(k_x+1,k_y+1);
        
        if k_c==k_diag+prom_base(seq_x(k_x),seq_y(k_y),score_m,score_mm) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,base_y(k_y)];
            k_x=k_x-1;
            k_y=k_y-1;
        elseif k_c==k_gx+gap_penalty(1) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,45];
            k_x=k_x-1;
        elseif k_c==k_gy+gap_penalty(1) then
            out_x = [out_x,45];
            out_y = [out_y,base_y(k_y)];
            k_y=k_y-1;
        elseif k_c==0;
            break;
        end
    end
    
    // Write the output
    if (k_x>0&k_y==0) then
        pre_x = ascii(base_x(1:k_x));
        pre_y = ascii(45*ones(1,k_x));
        align_x = strcat([pre_x,strrev(ascii(out_x)),post_x]);
        align_y = strcat([pre_y,strrev(ascii(out_y)),post_y]);
    elseif (k_x==0&k_y>0) then
        pre_y = ascii(base_y(1:k_y));
        pre_x = ascii(45*ones(1,k_y));
        align_x = strcat([pre_x,strrev(ascii(out_x)),post_x]);
        align_y = strcat([pre_y,strrev(ascii(out_y)),post_y]);       
    else
        align_x = strrev(ascii(out_x));
        align_y = strrev(ascii(out_y));          
    end
    prom_pos = length(pre_y);     
endfunction

function matrix_result = scorematrix_semiglobal(seq_x,seq_y,score_m,score_mm,gap_penalty)
    //Function to perform local search
    len_x = length(seq_x);
    len_y = length(seq_y);
    basic_mat = zeros(len_x+1,len_y+1);
    
    //No need to initialize because gap penalty will always be negative
    //Recurrance relation
    for k_x=2:(len_x+1)
        for k_y=2:(len_y+1)
            score_match = basic_mat(k_x-1,k_y-1)+comp_base(seq_x(k_x-1),seq_y(k_y-1),score_m,score_mm);
            score_gap_x = basic_mat(k_x,k_y-1)+gap_penalty(1);
            score_gap_y = basic_mat(k_x-1,k_y)+gap_penalty(1);
            basic_mat(k_x,k_y)=max([score_match score_gap_x score_gap_y]); // Same as global search
        end
    end
    matrix_result = basic_mat;
endfunction

function [align_x,align_y] = traceback_semiglobal(seq_x,seq_y,score_m,score_mm,gap_penalty)
    // Traceback for semi-global alignment
    
    // Assuming x is longer than y to make code simpler
    
    score_matrix = scorematrix_semiglobal(seq_x,seq_y,score_m,score_mm,gap_penalty);
        
    len_x = length(seq_x);
    len_y = length(seq_y);
    
    base_x = seq_x;
    base_y = seq_y;
    
    out_x = [];
    out_y = [];
    
    k_x = len_x;
    [n_max,k_y] = max(score_matrix(len_x+1,:)); 
    k_y=k_y-1;
    
    post_x = ascii(seq_x((k_x+1):len_x));
    post_y = ascii(seq_y((k_y+1):len_y));
    
    // Perform the traceback
    while (k_x>0&k_y>0) then
        k_gx = score_matrix(k_x,k_y+1);
        k_gy = score_matrix(k_x+1,k_y);
        k_diag = score_matrix(k_x,k_y);
        k_c = score_matrix(k_x+1,k_y+1);
        
        if k_c==k_diag+comp_base(seq_x(k_x),seq_y(k_y),score_m,score_mm) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,base_y(k_y)];
            k_x=k_x-1;
            k_y=k_y-1;
        elseif k_c==k_gx+gap_penalty(1) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,45];
            k_x=k_x-1;
        elseif k_c==k_gy+gap_penalty(1) then
            out_x = [out_x,45];
            out_y = [out_y,base_y(k_y)];
            k_y=k_y-1;
        end
    end
    
    // Write the output
    if (k_x>0&k_y==0) then
        pre_x = ascii(base_x(1:k_x));
        pre_y = ascii(45*ones(1,k_x));
        align_x = strcat([pre_x,strrev(ascii(out_x)),post_x]);
        align_y = strcat([pre_y,strrev(ascii(out_y)),post_y]);
    elseif (k_x==0&k_y>0) then
        pre_y = ascii(base_y(1:k_y));
        pre_x = ascii(45*ones(1,k_y));
        align_x = strcat([pre_x,strrev(ascii(out_x)),post_x]);
        align_y = strcat([pre_y,strrev(ascii(out_y)),post_y]);       
    else
        align_x = strrev(ascii(out_x));
        align_y = strrev(ascii(out_y));          
    end     
endfunction

function result = remove_eols(text_in)
    // Remove EOLs of fasta file
    keys = find(text_in==10);
    if (isempty(keys)) then
        result = text_in;
    else
    text_out = text_in(1:(keys(1)-1));
    k_n = length(keys)-1;
    for k=1:k_n
        text_i = text_in((keys(k)+1):(keys(k+1)-1));
        text_out = [text_out, text_i];
    end
    result = [text_out,text_in((keys(k_n+1)+1):length(text_in))];
    end
endfunction

function comp=get_comp(base)
    if (base==65) then
        comp = 84;
    elseif (base==67) then
        comp = 71;
    elseif (base==71) then
        comp = 67;
    elseif (base==84) then
        comp = 65;
    end
endfunction

function bk=get_base_key(base)
    // Base key as A=1, C=2, G=3, T=4 (alphabetic)
    if (base==65) then
        bk=1;
    elseif (base==67) then
        bk=2;
    elseif (base==71) then
        bk=3;
    elseif (base==84) then
        bk=4;
    // For multiple possibilities assign one for consistency
    // Evenly distributed as much as possible A=3,C=3, G=3, T=2
    elseif (base==ascii('R')) then
        bk=1; // Assign A (A or G)    
    elseif (base==ascii('Y')) then
        bk=2; // Assign C (C or T) 
    elseif (base==ascii('S')) then
        bk=3; // Assign G (C or G) 
    elseif (base==ascii('W')) then
        bk=4; // Assign T (A or T) 
    elseif (base==ascii('K')) then
        bk=3; // Assign G (G or T) 
    elseif (base==ascii('M')) then
        bk=1; // Assign A (A or C) 
    elseif (base==ascii('B')) then
        bk=2; // Assign C (C or G or T) 
    elseif (base==ascii('D')) then
        bk=3; // Assign G (A or G or T) 
    elseif (base==ascii('H')) then
        bk=4; // Assign T (A or C or T) 
    elseif (base==ascii('V')) then
        bk=1; // Assign A (A or C or G)
    elseif (base==ascii('N')) then
        bk=2; // Assign C (any base)
    end
endfunction

function base_hist=get_base_hist(seq)
    hist=ones(1,4)/4;
    l_seq = length(seq);
    for pos=1:l_seq
        key=get_base_key(seq(pos));
        hist(key)=hist(key)+1;
    end
    base_hist=hist;
endfunction

function gen_code=get_fasta_at(file_name,g_pos,g_end,strand)
    //Estimate the necessary overread to compensate for the EOL charactors of FASTA files
    g_len = g_end-g_pos;
    if g_len>0 then
           n_extra = floor(g_len/70);
           n_offset = floor(g_pos/70);
           file_details = fileinfo(file_name);
           file_len = file_details(1);
           fd = mopen(file_name,'rb');
           mseek(0,fd);
           header = mgetl(fd,1);
           g_start = length(header);
           mseek(g_start+g_pos+n_offset,fd);
           raw_code = mget(g_len+n_extra,'c',fd);
           mclose(fd);
           code_i = remove_eols(raw_code);
           if strand==1 then
               gen_code = code_i;
           else
               //get complementary strand
               len = length(code_i);
               code_c = [];
               for k=1:len
                   code_c = [code_c,get_comp(code_i(k))];
                   gen_code = code_c;
               end
           end
    else
        gen_code = [];
    end
 
endfunction

function [gene_array_p,gene_array_n,noncoding_array_p,noncoding_array_n]=get_protein_pos_array(filename)
    // Get the coding and non-coding DNA positions from a protein table
    fd = mopen(filename,'r');
    data = mgetl(fd,1);
    ga_p = [];
    ga_n = [];
    
    nca_p = [];
    nca_n = [];
    
    nc_prev_p = 0;
    nc_prev_n = 0;
    
    while (~meof(fd)) 
        data = mgetl(fd,1);
        //disp(type(data));
        keys = strindex(data,ascii(9));
        if (isempty(keys)) then
            break;
        end
        p_data = strsplit(data,keys);
        pg_start = strtod(p_data(3));
        pg_stop = strtod(p_data(4));
        //disp(strcmp(p_data(5),'-'));
        //disp(strcmp(p_data(5),'+'));
        if (~isempty(p_data(5))) then
            if (strcmp(p_data(5),'-')==1) then
                ga_n = [ga_n; pg_start,pg_stop];
                nca_n = [nca_n; nc_prev_n, (pg_start-1)];
                nc_prev_n = pg_stop+1;
            else             
                ga_p = [ga_p; pg_start,pg_stop];
                nca_p = [nca_p; nc_prev_p, (pg_start-1)];
                nc_prev_p = pg_stop+1;
            end
        end
    end
    mclose(fd);
    
    gene_array_p=ga_p;
    noncoding_array_p=nca_p;
    gene_array_n=ga_n;
    noncoding_array_n=nca_n;
endfunction

function pos_freq_mat_out = update_pos_freq_mat(seq,pos_freq_mat_in,set_stop)
    [m,n]=size(pos_freq_mat_in);
    s_len = length(seq);
    stop_point = min([n s_len set_stop]); // Where to stop
    for k=1:stop_point
        pos_base_key = get_base_key(seq(k));
        pos_freq_mat_in(pos_base_key,k)=pos_freq_mat_in(pos_base_key,k)+1;
    end
    pos_freq_mat_out = pos_freq_mat_in;
endfunction

function ppm=get_ppm(pos_freq_matrix)
    // Get the position probability matrix
    [m,n]=size(pos_freq_matrix);
    col_sum = sum(pos_freq_matrix,1);
    ppm_temp = [];
    for k=1:n
        ppm_temp = [ppm_temp,pos_freq_matrix(:,k)/col_sum(k)];
    end
    ppm = ppm_temp;
endfunction

function [w,su]=ppm_info(ppm,p0)
    [m,n]=size(ppm);
    p_mat = [];
    for k=1:m
        p_row = ones(1,n)/p0(k);
        p_mat = [p_mat; p_row];
    end
    w = log(ppm.*p_mat)/log(2);
    su = ppm.*w;
endfunction