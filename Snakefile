configfile: "config.yaml"
Genomes = config["Genomes"]

indx_suff ='invalid'
if config["Query_Type"] == "search":
    indx_suff = 'rns' # rank not supported
elif config["Query_Type"] == "rank":
    indx_suff = 'rs' #rank supported

rule all:
    input:
        # dictionary uniform knot encoding
        expand("{fn}/{fn}.index_"+config["Indx_Type"]+"_"+config["Eps"]+"_"+
            config["Kmer_Size"]+"_"+indx_suff+".bin", fn=Genomes),
        expand("{fn}/{fn}.query_time_"+config["Indx_Type"]+"_"+config["Eps"]+"_"
            +config["Kmer_Size"]+"_"+indx_suff+".txt", fn=Genomes)

rule build_dictionary:
    input:
        genomeF = "{fn}/{fn}.ConcatenatedGenome.txt", #only bases concatenated
        saF = "{fn}/{fn}.sa.bin"
    params:
        kmer_size = config["Kmer_Size"],
        codePath = config["CodePath"],
        execPath = config["ExecPath"],
        flags = config["Compile_Flag"],
        lp_bits = config["LP_Bits"],
        epsilon = config["Eps"],
        indx_type = config["Indx_Type"],
        query_type = config["Query_Type"],
        sdsl_include_path = config["SDSL_Inc"],
        sdsl_lib_path = config["SDSL_Lib"],
    output:
        index_fn = "{fn}/{fn}.index_"+config["Indx_Type"]+"_"+config["Eps"]+"_"+
            config["Kmer_Size"]+"_"+indx_suff+".bin"
    shell:
        "g++ -std=c++17 {params.flags} -I {params.sdsl_include_path} -L {params.sdsl_lib_path} {params.codePath}/build_index.cpp "
        "{params.codePath}/BitPacking.cpp {params.codePath}/pla_index.cpp "
        "-o {params.execPath}/build_index "
        "-lsdsl -ldivsufsort -ldivsufsort64;"
        "/usr/bin/time -f \"%M,%e,%U,%S\" --output-file=memkb_sec_Usec_Ksec_dict_build.txt "
        "{params.execPath}/build_index {input.genomeF} {input.saF} "
        "{params.kmer_size} {params.epsilon} "
        "{output.index_fn} {params.lp_bits} "
        "{params.indx_type} {params.query_type} "

rule benchmark:
    input:
        index_fn = "{fn}/{fn}.index_"+config["Indx_Type"]+"_"+config["Eps"]+"_"+
            config["Kmer_Size"]+"_"+indx_suff+".bin",
        saF = "{fn}/{fn}.sa.bin",
        query_file = "{fn}/{fn}."+config["Query_Tag"]+".query.txt",
        genomeF = "{fn}/{fn}.ConcatenatedGenome.txt" #only bases concatenated
    params:
        knot_bs_thres = config["Knot_BS_Thres"],
        kmer_size = config["Kmer_Size"],
        codePath = config["CodePath"],
        execPath = config["ExecPath"],
        flags = config["Compile_Flag"],
        indx_type = config["Indx_Type"],
        query_type = config["Query_Type"],
        sdsl_include_path = config["SDSL_Inc"],
        sdsl_lib_path = config["SDSL_Lib"]

    output:
        runInfo = "{fn}/{fn}.query_time_"+config["Indx_Type"]+"_"+config["Eps"]+"_"
            +config["Kmer_Size"]+"_"+indx_suff+".txt"
    shell:
        "g++ -std=c++17 {params.flags} -I {params.sdsl_include_path} -L {params.sdsl_lib_path} {params.codePath}/benchmark_index.cpp "
        "{params.codePath}/BitPacking.cpp {params.codePath}/pla_index.cpp "
        "-o {params.execPath}/benchmark_index "
        "-lsdsl -ldivsufsort -ldivsufsort64;"
        "/usr/bin/time -f \"%M,%e,%U,%S\" --output-file=memkb_sec_Usec_Ksec_query.txt "
        "{params.execPath}/benchmark_index {input.saF} {input.genomeF} {params.kmer_size} {input.query_file} "
        "{input.index_fn}  {output.runInfo} {params.knot_bs_thres} "
        "{params.indx_type} {params.query_type} "
