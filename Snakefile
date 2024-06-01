configfile: "config.yaml"
Genomes = config["Genomes"]
Genomes = [gn.rstrip('/') for gn in Genomes]

fast_rank_flag = " "
if config["Fast_Rank"] != 'n':
    fast_rank_flag = " -r "

rule all:
    input:
        # dictionary uniform knot encoding
        expand("{fn}/{fn}.index", fn=Genomes),
        expand("{fn}/{fn}.query_time.txt", fn=Genomes)

rule build:
    input:
        genomeF = "{fn}/{fn}.processed.fasta", #single header line, bases concatenated ideally
        saF = "{fn}/{fn}.sa.bin"
    params:
        kmer_size = config["Kmer_Size"],
        codePath = config["CodePath"],
        execPath = config["ExecPath"],
        flags = config["Compile_Flag"],
        lookup = config["Lookup"],
        epsilon = config["Eps"],
        indx_type = config["Indx_Type"],
        is_fast_rank_enabled = fast_rank_flag,
        sdsl_include_path = config["SDSL_Inc"],
        sdsl_lib_path = config["SDSL_Lib"],
    output:
        index_fn = "{fn}/{fn}.index"
    shell:
        "g++ -std=c++17 {params.flags} -I {params.sdsl_include_path} -L {params.sdsl_lib_path} {params.codePath}/build_pla_index.cpp "
        "{params.codePath}/pla_index.cpp "
        "{params.codePath}/cmdline.cpp "
        "-o {params.execPath}/build_pla_index "
        "-lsdsl -ldivsufsort -ldivsufsort64;"
        # "/usr/bin/time -f \"%M,%e,%U,%S\" --output-file=memkb_sec_Usec_Ksec_dict_build.txt "
        "{params.execPath}/build_pla_index -g {input.genomeF} -s {input.saF} "
        "-k {params.kmer_size} -e {params.epsilon} "
        "-o {output.index_fn} -l {params.lookup} "
        "-t {params.indx_type} {params.is_fast_rank_enabled} "

rule query:
    input:
        index_fn = "{fn}/{fn}.index",
        saF = "{fn}/{fn}.sa.bin",
        query_file = "{fn}/{fn}."+config["Query_Tag"]+".query.txt",
        genomeF = "{fn}/{fn}.processed.fasta"
    params:
        codePath = config["CodePath"],
        execPath = config["ExecPath"],
        flags = config["Compile_Flag"],
        query_type = config["Query_Type"],
        sdsl_include_path = config["SDSL_Inc"],
        sdsl_lib_path = config["SDSL_Lib"]
    output:
        runInfo = "{fn}/{fn}.query_time.txt"
    shell:
        "g++ -std=c++17 {params.flags} -I {params.sdsl_include_path} -L {params.sdsl_lib_path} {params.codePath}/query_pla_index.cpp "
        " {params.codePath}/pla_index.cpp "
        "{params.codePath}/cmdline.cpp "
        "-o {params.execPath}/query_pla_index "
        "-lsdsl -ldivsufsort -ldivsufsort64;"
        # "/usr/bin/time -f \"%M,%e,%U,%S\" --output-file=memkb_sec_Usec_Ksec_query.txt "
        "{params.execPath}/query_pla_index -g {input.genomeF} -s {input.saF} -q {input.query_file} "
        "-i {input.index_fn} --run_info {output.runInfo} "
        "--query_type {params.query_type} "
