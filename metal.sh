# ===============================================================
# Describe: calculated META by metal
# Author: Weber
# Date: 2025.08.13
# Update: 2025.08.13
# version: 1.0
# Parameters: typing "bash metal.sh --help"
# ===============================================================
#! /usr/bin/bash


Help()
{
	# Display Help
	echo
	echo "***Pipeline for metal***"
	echo "Add description of the script functions here:"
	echo
    echo "  -p | --pheno       : Phenotype name, used for output file naming."
	echo "  -m | --mvp         : Summary statisitcs file of mvp_meta."
	echo "  -g | --gwas        : Summary statisitcs file of GWAS."
    echo "  -t | --target      : Target sample prefix for PRS-CS after GWAS (Optional)."
	echo "  -v | --valid       : Validation sample prefix for PRS-CS after GWAS (Optional)."
	echo 
}


metafor_path='/home/Weber/Metal/metafor.r'
metal_r_path='/home/Weber/Metal/gwas_metal.r'
meta_merge_path='/home/Weber/Pipeline/META/meta_merge.py'


re='^(--help|-h)$'
if [[ $1 =~ $re ]];
then
	Help
else
	while [ "$#" -gt 0 ]; do
		case "$1" in
            -p|--pheno) pheno="$2"; shift 2;;
            -m|--mvp) mvp="$2"; shift 2;;
			-g|--gwas) gwas="$2"; shift 2;;
            -t|--target) target="$2"; shift 2;;
            -v|--valid) valid="$2"; shift 2;;
			*) echo "unknown option: $1" >&2; exit 1;;
		esac
	done

    echo -e "Activating conda environment: myenvR"
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate myenvR

    eval Rscript $metal_r_path $mvp $gwas
    # eval Rscript $metafor_path
    echo -e "Metal analysis completed."


    echo -e "Plotting Manhattan plot using RunTopR_meta.R..."
    
    
    Rscript /home/Weber/Pipeline/META/RunTopR_meta.R \
        --input metal_input.txt \
        --meta meta_fixed_results.txt \
        --output ${pheno}_metal_GWAS_results_fixed \
        --thresh 1e-5 \
        --annotate 5e-8 \
        --ymax NA 

    Rscript /home/Weber/Pipeline/META/RunTopR_meta.R \
        --input metal_input.txt \
        --meta meta_random_results.txt \
        --output ${pheno}_metal_GWAS_results_random \
        --thresh 1e-5 \
        --annotate 5e-8 \
        --ymax NA 


    echo -e "Plot manhattan had done..."

    conda activate myenv

    echo -e "Running PRS-CSX for target and validation sets..."
    # 若指定 --target 則準備執行 PRSCS
    if [ ! -z "$target" ] && [ ! -z "$valid" ]; then
        echo -e "[INFO] GWAS done. Will run PRSCS using target: $target"
        echo -e "[INFO] Valid set: $valid"

        # 需調整
        prscsx_cmd="bash /home/Weber/Pipeline/PRS/PRSCSX_for_target.sh \
                    -b /SNParray/SourceShare/20240321_50w_Imputation/step12-pgen2bed/Axiom_imputed_r2.MAF \
                    -g finaloutput_fixed_input.txt \
                    -v v2 \
                    -n 500000 \
                    -c all \
                    -o meta_prs/${pheno}_fixed_META_PRS \
                    --target_list $target \
                    --val_list $valid"
                    #--pop_list $pop_list"
                    
        echo -e "[INFO] Running PRSCSX with the following command:"
        echo "$prscsx_cmd"
        eval $prscsx_cmd

        # echo -e "[INFO] finish fixed_meta PRSCSX"
        #  # 需調整
        # prscsx_cmd="bash /home/Weber/Pipeline/PRS/PRSCSX_for_target.sh \
        #             -b /SNParray/SourceShare/20240321_50w_Imputation/step12-pgen2bed/Axiom_imputed_r2.MAF \
        #             -g finaloutput_random_input.txt \
        #             -v v2 \
        #             -n 500000 \
        #             -c all \
        #             -o meta_prs/${pheno}_random_META_PRS \
        #             --target_list $target \
        #             --val_list $valid"
                    # --pop_list $pop_list"

                    
        # echo -e "[INFO] Running PRSCSX with the following command:"
        # echo "$prscsx_cmd"
        # eval $prscsx_cmd
        
    fi
    echo -e "PRS-CSX analysis completed."

    echo -e "Merge fixed and random meta..."
    eval python $meta_merge_path "${pheno}"
    echo -e "Merge meta completed."

fi