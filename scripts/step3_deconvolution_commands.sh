#!/usr/bin/env bash
# CTSE benchmark command catalog: step 3 deconvolution
#
# Run selected commands from scripts/. This file is a reproducibility record,
# not a single end-to-end job: commands can refresh existing benchmark outputs.
# BLUE and scTAPE require separately distributed processed h5ad inputs.
# CIBERSORTx requires credentials and a licensed container supplied at run time.
# Set CIBERSORTX_USERNAME, CIBERSORTX_TOKEN, and CIBERSORTX_CONTAINER_PATH
# before selecting a CIBERSORTx command; secrets are not stored here.

# =============================================================================
# NOISY BULK INPUT PREPARATION
# Regenerate the configured noisy bulk inputs and their QC summaries (seed 123).
# =============================================================================
# BRCA_Bassez2021 — noisy bulk inputs (optional)
Rscript step4_run_add_bulk_noise.R --dataset BRCA_Bassez2021 --seed 123
# CRC_Pelka2021 — noisy bulk inputs
Rscript step4_run_add_bulk_noise.R --dataset CRC_Pelka2021 --seed 123
# LUAD_Kim2020 — noisy bulk inputs
Rscript step4_run_add_bulk_noise.R --dataset LUAD_Kim2020 --seed 123
# PBMC_Perez2022 — noisy bulk inputs
Rscript step4_run_add_bulk_noise.R --dataset PBMC_Perez2022 --seed 123
# ROSMAP_AD92_Xiong2023 — noisy bulk inputs
Rscript step4_run_add_bulk_noise.R --dataset ROSMAP_AD92_Xiong2023 --seed 123

# Deconvolution commands below reproduce the retained deconv_res method folders.

# =============================================================================
# DECONVOLUTION — BRCA_Bassez2021
# Configurations: config01–config09, plus config02_tmm and config02_uq.
# =============================================================================
Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset BRCA_Bassez2021 --config_id config01 --n_core 15
python3 ../DALE_Eval/run/run_BLUE.py --dataset BRCA_Bassez2021 --config_id config01 --n_core 15
python3 ../DALE_Eval/run/run_scTAPE.py --dataset BRCA_Bassez2021 --config_id config01 --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset BRCA_Bassez2021 --config_id config01 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset BRCA_Bassez2021 --config_id config01 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset BRCA_Bassez2021 --config_id config01 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset BRCA_Bassez2021 --config_id config01 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset BRCA_Bassez2021 --config_id config01 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset BRCA_Bassez2021 --config_id config01 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset BRCA_Bassez2021 --config_id config01 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset BRCA_Bassez2021 --config_id config02 --n_core 15
python3 ../DALE_Eval/run/run_BLUE.py --dataset BRCA_Bassez2021 --config_id config02 --n_core 15
python3 ../DALE_Eval/run/run_scTAPE.py --dataset BRCA_Bassez2021 --config_id config02 --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset BRCA_Bassez2021 --config_id config02 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset BRCA_Bassez2021 --config_id config02 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset BRCA_Bassez2021 --config_id config02 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset BRCA_Bassez2021 --config_id config02 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset BRCA_Bassez2021 --config_id config02 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset BRCA_Bassez2021 --config_id config02 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset BRCA_Bassez2021 --config_id config02 --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset BRCA_Bassez2021 --config_id config02_tmm --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"

Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset BRCA_Bassez2021 --config_id config02_tmm --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset BRCA_Bassez2021 --config_id config02_tmm --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset BRCA_Bassez2021 --config_id config02_tmm --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset BRCA_Bassez2021 --config_id config02_tmm --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset BRCA_Bassez2021 --config_id config02_tmm --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset BRCA_Bassez2021 --config_id config02_tmm --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset BRCA_Bassez2021 --config_id config02_uq --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"

Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset BRCA_Bassez2021 --config_id config02_uq --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset BRCA_Bassez2021 --config_id config02_uq --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset BRCA_Bassez2021 --config_id config02_uq --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset BRCA_Bassez2021 --config_id config02_uq --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset BRCA_Bassez2021 --config_id config02_uq --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset BRCA_Bassez2021 --config_id config02_uq --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset BRCA_Bassez2021 --config_id config03 --n_core 15
python3 ../DALE_Eval/run/run_BLUE.py --dataset BRCA_Bassez2021 --config_id config03 --n_core 15
python3 ../DALE_Eval/run/run_scTAPE.py --dataset BRCA_Bassez2021 --config_id config03 --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset BRCA_Bassez2021 --config_id config03 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset BRCA_Bassez2021 --config_id config03 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset BRCA_Bassez2021 --config_id config03 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset BRCA_Bassez2021 --config_id config03 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset BRCA_Bassez2021 --config_id config03 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset BRCA_Bassez2021 --config_id config03 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset BRCA_Bassez2021 --config_id config03 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset BRCA_Bassez2021 --config_id config04 --n_core 15
python3 ../DALE_Eval/run/run_BLUE.py --dataset BRCA_Bassez2021 --config_id config04 --n_core 15
python3 ../DALE_Eval/run/run_scTAPE.py --dataset BRCA_Bassez2021 --config_id config04 --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset BRCA_Bassez2021 --config_id config04 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset BRCA_Bassez2021 --config_id config04 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset BRCA_Bassez2021 --config_id config04 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset BRCA_Bassez2021 --config_id config04 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset BRCA_Bassez2021 --config_id config04 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset BRCA_Bassez2021 --config_id config04 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset BRCA_Bassez2021 --config_id config04 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset BRCA_Bassez2021 --config_id config05 --n_core 15
python3 ../DALE_Eval/run/run_BLUE.py --dataset BRCA_Bassez2021 --config_id config05 --n_core 15
python3 ../DALE_Eval/run/run_scTAPE.py --dataset BRCA_Bassez2021 --config_id config05 --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset BRCA_Bassez2021 --config_id config05 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset BRCA_Bassez2021 --config_id config05 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset BRCA_Bassez2021 --config_id config05 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset BRCA_Bassez2021 --config_id config05 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset BRCA_Bassez2021 --config_id config05 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset BRCA_Bassez2021 --config_id config05 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset BRCA_Bassez2021 --config_id config05 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset BRCA_Bassez2021 --config_id config06 --n_core 15
python3 ../DALE_Eval/run/run_BLUE.py --dataset BRCA_Bassez2021 --config_id config06 --n_core 15
python3 ../DALE_Eval/run/run_scTAPE.py --dataset BRCA_Bassez2021 --config_id config06 --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset BRCA_Bassez2021 --config_id config06 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset BRCA_Bassez2021 --config_id config06 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset BRCA_Bassez2021 --config_id config06 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset BRCA_Bassez2021 --config_id config06 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset BRCA_Bassez2021 --config_id config06 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset BRCA_Bassez2021 --config_id config06 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset BRCA_Bassez2021 --config_id config06 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset BRCA_Bassez2021 --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset BRCA_Bassez2021 --config_id config07 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset BRCA_Bassez2021 --config_id config07 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset BRCA_Bassez2021 --config_id config07 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset BRCA_Bassez2021 --config_id config07 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset BRCA_Bassez2021 --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset BRCA_Bassez2021 --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset BRCA_Bassez2021 --config_id config07 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset BRCA_Bassez2021 --config_id config08 --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset BRCA_Bassez2021 --config_id config08 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset BRCA_Bassez2021 --config_id config08 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset BRCA_Bassez2021 --config_id config08 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset BRCA_Bassez2021 --config_id config08 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset BRCA_Bassez2021 --config_id config08 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset BRCA_Bassez2021 --config_id config08 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset BRCA_Bassez2021 --config_id config08 --n_core 15

Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset BRCA_Bassez2021 --config_id config09 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset BRCA_Bassez2021 --config_id config09 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset BRCA_Bassez2021 --config_id config09 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset BRCA_Bassez2021 --config_id config09 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset BRCA_Bassez2021 --config_id config09 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset BRCA_Bassez2021 --config_id config09 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset BRCA_Bassez2021 --config_id config09 --n_core 15



# =============================================================================
# DECONVOLUTION — CRC_Pelka2021
# Configurations: config01, config02, config07, config08, and config09.
# =============================================================================
Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset CRC_Pelka2021 --config_id config01 --n_core 15
python3 ../DALE_Eval/run/run_BLUE.py --dataset CRC_Pelka2021 --config_id config01 --n_core 15
python3 ../DALE_Eval/run/run_scTAPE.py --dataset CRC_Pelka2021 --config_id config01 --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset CRC_Pelka2021 --config_id config01 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset CRC_Pelka2021 --config_id config01 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset CRC_Pelka2021 --config_id config01 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset CRC_Pelka2021 --config_id config01 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset CRC_Pelka2021 --config_id config01 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset CRC_Pelka2021 --config_id config01 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset CRC_Pelka2021 --config_id config01 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset CRC_Pelka2021 --config_id config02 --n_core 15
python3 ../DALE_Eval/run/run_BLUE.py --dataset CRC_Pelka2021 --config_id config02 --n_core 15
python3 ../DALE_Eval/run/run_scTAPE.py --dataset CRC_Pelka2021 --config_id config02 --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset CRC_Pelka2021 --config_id config02 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset CRC_Pelka2021 --config_id config02 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset CRC_Pelka2021 --config_id config02 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset CRC_Pelka2021 --config_id config02 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset CRC_Pelka2021 --config_id config02 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset CRC_Pelka2021 --config_id config02 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset CRC_Pelka2021 --config_id config02 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset CRC_Pelka2021 --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset CRC_Pelka2021 --config_id config07 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset CRC_Pelka2021 --config_id config07 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset CRC_Pelka2021 --config_id config07 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset CRC_Pelka2021 --config_id config07 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset CRC_Pelka2021 --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset CRC_Pelka2021 --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset CRC_Pelka2021 --config_id config07 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset CRC_Pelka2021 --config_id config08 --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset CRC_Pelka2021 --config_id config08 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset CRC_Pelka2021 --config_id config08 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset CRC_Pelka2021 --config_id config08 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset CRC_Pelka2021 --config_id config08 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset CRC_Pelka2021 --config_id config08 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset CRC_Pelka2021 --config_id config08 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset CRC_Pelka2021 --config_id config08 --n_core 15

Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset CRC_Pelka2021 --config_id config09 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset CRC_Pelka2021 --config_id config09 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset CRC_Pelka2021 --config_id config09 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset CRC_Pelka2021 --config_id config09 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset CRC_Pelka2021 --config_id config09 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset CRC_Pelka2021 --config_id config09 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset CRC_Pelka2021 --config_id config09 --n_core 15



# =============================================================================
# DECONVOLUTION — LUAD_Kim2020
# Configurations: config01, config02, config07, config08, and config09.
# =============================================================================
Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset LUAD_Kim2020 --config_id config01 --n_core 15
python3 ../DALE_Eval/run/run_BLUE.py --dataset LUAD_Kim2020 --config_id config01 --n_core 15 --extra_args 'samplenum_per_ct=10;val_samplenum_per_patient=4'
python3 ../DALE_Eval/run/run_scTAPE.py --dataset LUAD_Kim2020 --config_id config01 --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset LUAD_Kim2020 --config_id config01 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset LUAD_Kim2020 --config_id config01 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset LUAD_Kim2020 --config_id config01 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset LUAD_Kim2020 --config_id config01 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset LUAD_Kim2020 --config_id config01 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset LUAD_Kim2020 --config_id config01 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset LUAD_Kim2020 --config_id config01 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset LUAD_Kim2020 --config_id config02 --n_core 15
python3 ../DALE_Eval/run/run_BLUE.py --dataset LUAD_Kim2020 --config_id config02 --n_core 15 --extra_args 'samplenum_per_ct=10;val_samplenum_per_patient=4'
python3 ../DALE_Eval/run/run_scTAPE.py --dataset LUAD_Kim2020 --config_id config02 --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset LUAD_Kim2020 --config_id config02 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset LUAD_Kim2020 --config_id config02 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset LUAD_Kim2020 --config_id config02 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset LUAD_Kim2020 --config_id config02 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset LUAD_Kim2020 --config_id config02 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset LUAD_Kim2020 --config_id config02 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset LUAD_Kim2020 --config_id config02 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset LUAD_Kim2020 --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset LUAD_Kim2020 --config_id config07 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset LUAD_Kim2020 --config_id config07 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset LUAD_Kim2020 --config_id config07 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset LUAD_Kim2020 --config_id config07 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset LUAD_Kim2020 --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset LUAD_Kim2020 --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset LUAD_Kim2020 --config_id config07 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset LUAD_Kim2020 --config_id config08 --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset LUAD_Kim2020 --config_id config08 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset LUAD_Kim2020 --config_id config08 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset LUAD_Kim2020 --config_id config08 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset LUAD_Kim2020 --config_id config08 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset LUAD_Kim2020 --config_id config08 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset LUAD_Kim2020 --config_id config08 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset LUAD_Kim2020 --config_id config08 --n_core 15

Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset LUAD_Kim2020 --config_id config09 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset LUAD_Kim2020 --config_id config09 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset LUAD_Kim2020 --config_id config09 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset LUAD_Kim2020 --config_id config09 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset LUAD_Kim2020 --config_id config09 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset LUAD_Kim2020 --config_id config09 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset LUAD_Kim2020 --config_id config09 --n_core 15



# =============================================================================
# DECONVOLUTION — PBMC_1k1k
# Configurations: config01 and config07.
# =============================================================================
Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset PBMC_1k1k --config_id config01 --n_core 15
python3 ../DALE_Eval/run/run_BLUE.py --dataset PBMC_1k1k --config_id config01 --n_core 15 --extra_args 'ref_sample_subset=true;ref_sample_subset_n=40'
python3 ../DALE_Eval/run/run_scTAPE.py --dataset PBMC_1k1k --config_id config01 --n_core 15 --extra_args 'ref_sample_subset=true;ref_sample_subset_n=40'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_1k1k --config_id config01 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_1k1k --config_id config01 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset PBMC_1k1k --config_id config01 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset PBMC_1k1k --config_id config01 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset PBMC_1k1k --config_id config01 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset PBMC_1k1k --config_id config01 --n_core 15 --extra_args 'use_limma_top_genes=true;top_n=3000;chunk_size=300'

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset PBMC_1k1k --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_1k1k --config_id config07 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_1k1k --config_id config07 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset PBMC_1k1k --config_id config07 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset PBMC_1k1k --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset PBMC_1k1k --config_id config07 --n_core 15



# =============================================================================
# DECONVOLUTION — PBMC_Perez2022
# Configurations: config01–config09.
# =============================================================================
Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset PBMC_Perez2022 --config_id config01 --n_core 15
python3 ../DALE_Eval/run/run_BLUE.py --dataset PBMC_Perez2022 --config_id config01 --n_core 15 --extra_args 'ref_sample_subset=true;ref_sample_subset_n=40'
python3 ../DALE_Eval/run/run_scTAPE.py --dataset PBMC_Perez2022 --config_id config01 --n_core 15 --extra_args 'ref_sample_subset=true;ref_sample_subset_n=40'
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset PBMC_Perez2022 --config_id config01 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_Perez2022 --config_id config01 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_Perez2022 --config_id config01 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset PBMC_Perez2022 --config_id config01 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset PBMC_Perez2022 --config_id config01 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset PBMC_Perez2022 --config_id config01 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset PBMC_Perez2022 --config_id config01 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset PBMC_Perez2022 --config_id config02 --n_core 15
python3 ../DALE_Eval/run/run_BLUE.py --dataset PBMC_Perez2022 --config_id config02 --n_core 15 --extra_args 'ref_sample_subset=true;ref_sample_subset_n=40'
python3 ../DALE_Eval/run/run_scTAPE.py --dataset PBMC_Perez2022 --config_id config02 --n_core 15 --extra_args 'ref_sample_subset=true;ref_sample_subset_n=40'
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset PBMC_Perez2022 --config_id config02 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_Perez2022 --config_id config02 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_Perez2022 --config_id config02 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset PBMC_Perez2022 --config_id config02 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset PBMC_Perez2022 --config_id config02 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset PBMC_Perez2022 --config_id config02 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset PBMC_Perez2022 --config_id config02 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset PBMC_Perez2022 --config_id config03 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset PBMC_Perez2022 --config_id config04 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset PBMC_Perez2022 --config_id config05 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset PBMC_Perez2022 --config_id config06 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset PBMC_Perez2022 --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset PBMC_Perez2022 --config_id config07 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_Perez2022 --config_id config07 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_Perez2022 --config_id config07 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset PBMC_Perez2022 --config_id config07 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset PBMC_Perez2022 --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset PBMC_Perez2022 --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset PBMC_Perez2022 --config_id config07 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset PBMC_Perez2022 --config_id config08 --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset PBMC_Perez2022 --config_id config08 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_Perez2022 --config_id config08 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_Perez2022 --config_id config08 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset PBMC_Perez2022 --config_id config08 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset PBMC_Perez2022 --config_id config08 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset PBMC_Perez2022 --config_id config08 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset PBMC_Perez2022 --config_id config08 --n_core 15

Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset PBMC_Perez2022 --config_id config09 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_Perez2022 --config_id config09 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_Perez2022 --config_id config09 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset PBMC_Perez2022 --config_id config09 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset PBMC_Perez2022 --config_id config09 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset PBMC_Perez2022 --config_id config09 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset PBMC_Perez2022 --config_id config09 --n_core 15



# =============================================================================
# DECONVOLUTION — PBMC_refined_1k1k
# Configurations: config01 and config07.
# =============================================================================
Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset PBMC_refined_1k1k --config_id config01 --n_core 15
python3 ../DALE_Eval/run/run_BLUE.py --dataset PBMC_refined_1k1k --config_id config01 --n_core 15 --extra_args 'ref_sample_subset=true;ref_sample_subset_n=40;cell_type_col=cell_type_refined'
python3 ../DALE_Eval/run/run_scTAPE.py --dataset PBMC_refined_1k1k --config_id config01 --n_core 15 --extra_args 'ref_sample_subset=true;ref_sample_subset_n=40;cell_type_col=cell_type_refined'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_refined_1k1k --config_id config01 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_refined_1k1k --config_id config01 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset PBMC_refined_1k1k --config_id config01 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false;use_limma_top_genes=true;top_n=3000'
Rscript ../DALE_Eval/run/run_TCA.R --dataset PBMC_refined_1k1k --config_id config01 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset PBMC_refined_1k1k --config_id config01 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset PBMC_refined_1k1k --config_id config01 --n_core 15 --extra_args 'use_limma_top_genes=true;top_n=3000;chunk_size=300'

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset PBMC_refined_1k1k --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_refined_1k1k --config_id config07 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_refined_1k1k --config_id config07 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_TCA.R --dataset PBMC_refined_1k1k --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset PBMC_refined_1k1k --config_id config07 --n_core 15



# =============================================================================
# DECONVOLUTION — PBMC_refined_Perez2022
# Configurations: config01, config07, and config09.
# =============================================================================
Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset PBMC_refined_Perez2022 --config_id config01 --n_core 15
python3 ../DALE_Eval/run/run_BLUE.py --dataset PBMC_refined_Perez2022 --config_id config01 --n_core 15 --extra_args 'ref_sample_subset=true;ref_sample_subset_n=40;cell_type_col=cell_type_refined'
python3 ../DALE_Eval/run/run_scTAPE.py --dataset PBMC_refined_Perez2022 --config_id config01 --n_core 15 --extra_args 'ref_sample_subset=true;ref_sample_subset_n=40;cell_type_col=cell_type_refined'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_refined_Perez2022 --config_id config01 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_refined_Perez2022 --config_id config01 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset PBMC_refined_Perez2022 --config_id config01 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false;use_limma_top_genes=true;top_n=3000'
Rscript ../DALE_Eval/run/run_TCA.R --dataset PBMC_refined_Perez2022 --config_id config01 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset PBMC_refined_Perez2022 --config_id config01 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset PBMC_refined_Perez2022 --config_id config01 --n_core 15 --extra_args 'use_limma_top_genes=true;top_n=3000;chunk_size=300'

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset PBMC_refined_Perez2022 --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_refined_Perez2022 --config_id config07 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_refined_Perez2022 --config_id config07 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset PBMC_refined_Perez2022 --config_id config07 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset PBMC_refined_Perez2022 --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset PBMC_refined_Perez2022 --config_id config07 --n_core 15

Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_refined_Perez2022 --config_id config09 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset PBMC_refined_Perez2022 --config_id config09 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset PBMC_refined_Perez2022 --config_id config09 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false;use_limma_top_genes=true;top_n=3000'
Rscript ../DALE_Eval/run/run_TCA.R --dataset PBMC_refined_Perez2022 --config_id config09 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset PBMC_refined_Perez2022 --config_id config09 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset PBMC_refined_Perez2022 --config_id config09 --n_core 15



# =============================================================================
# DECONVOLUTION — ROSMAP_AD430_Mathys2023
# Configurations: config01, config02, config07, and config08.
# =============================================================================
Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset ROSMAP_AD430_Mathys2023 --config_id config01 --n_core 15
python3 ../DALE_Eval/run/run_BLUE.py --dataset ROSMAP_AD430_Mathys2023 --config_id config01 --n_core 15 --extra_args 'ref_sample_subset=true;ref_sample_subset_n=25'
python3 ../DALE_Eval/run/run_scTAPE.py --dataset ROSMAP_AD430_Mathys2023 --config_id config01 --n_core 15 --extra_args 'ref_sample_subset=true;ref_sample_subset_n=25'
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset ROSMAP_AD430_Mathys2023 --config_id config01 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset ROSMAP_AD430_Mathys2023 --config_id config01 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset ROSMAP_AD430_Mathys2023 --config_id config01 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset ROSMAP_AD430_Mathys2023 --config_id config01 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset ROSMAP_AD430_Mathys2023 --config_id config01 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset ROSMAP_AD430_Mathys2023 --config_id config01 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset ROSMAP_AD430_Mathys2023 --config_id config01 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset ROSMAP_AD430_Mathys2023 --config_id config02 --n_core 15
python3 ../DALE_Eval/run/run_BLUE.py --dataset ROSMAP_AD430_Mathys2023 --config_id config02 --n_core 15 --extra_args 'ref_sample_subset=true;ref_sample_subset_n=25'
python3 ../DALE_Eval/run/run_scTAPE.py --dataset ROSMAP_AD430_Mathys2023 --config_id config02 --n_core 15 --extra_args 'ref_sample_subset=true;ref_sample_subset_n=25'
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset ROSMAP_AD430_Mathys2023 --config_id config02 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset ROSMAP_AD430_Mathys2023 --config_id config02 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset ROSMAP_AD430_Mathys2023 --config_id config02 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset ROSMAP_AD430_Mathys2023 --config_id config02 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset ROSMAP_AD430_Mathys2023 --config_id config02 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset ROSMAP_AD430_Mathys2023 --config_id config02 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset ROSMAP_AD430_Mathys2023 --config_id config02 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset ROSMAP_AD430_Mathys2023 --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset ROSMAP_AD430_Mathys2023 --config_id config07 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset ROSMAP_AD430_Mathys2023 --config_id config07 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset ROSMAP_AD430_Mathys2023 --config_id config07 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset ROSMAP_AD430_Mathys2023 --config_id config07 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset ROSMAP_AD430_Mathys2023 --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset ROSMAP_AD430_Mathys2023 --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset ROSMAP_AD430_Mathys2023 --config_id config07 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset ROSMAP_AD430_Mathys2023 --config_id config08 --n_core 15



# =============================================================================
# DECONVOLUTION — ROSMAP_AD92_Xiong2023
# Configurations: config01–config09.
# =============================================================================
Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset ROSMAP_AD92_Xiong2023 --config_id config01 --n_core 15
python3 ../DALE_Eval/run/run_BLUE.py --dataset ROSMAP_AD92_Xiong2023 --config_id config01 --n_core 15 --extra_args 'ref_sample_subset=true;ref_sample_subset_n=25'
python3 ../DALE_Eval/run/run_scTAPE.py --dataset ROSMAP_AD92_Xiong2023 --config_id config01 --n_core 15 --extra_args 'ref_sample_subset=true;ref_sample_subset_n=25'
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset ROSMAP_AD92_Xiong2023 --config_id config01 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset ROSMAP_AD92_Xiong2023 --config_id config01 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset ROSMAP_AD92_Xiong2023 --config_id config01 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset ROSMAP_AD92_Xiong2023 --config_id config01 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset ROSMAP_AD92_Xiong2023 --config_id config01 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset ROSMAP_AD92_Xiong2023 --config_id config01 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset ROSMAP_AD92_Xiong2023 --config_id config01 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset ROSMAP_AD92_Xiong2023 --config_id config02 --n_core 15
python3 ../DALE_Eval/run/run_BLUE.py --dataset ROSMAP_AD92_Xiong2023 --config_id config02 --n_core 15 --extra_args 'ref_sample_subset=true;ref_sample_subset_n=25'
python3 ../DALE_Eval/run/run_scTAPE.py --dataset ROSMAP_AD92_Xiong2023 --config_id config02 --n_core 15 --extra_args 'ref_sample_subset=true;ref_sample_subset_n=25'
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset ROSMAP_AD92_Xiong2023 --config_id config02 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset ROSMAP_AD92_Xiong2023 --config_id config02 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset ROSMAP_AD92_Xiong2023 --config_id config02 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset ROSMAP_AD92_Xiong2023 --config_id config02 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset ROSMAP_AD92_Xiong2023 --config_id config02 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset ROSMAP_AD92_Xiong2023 --config_id config02 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset ROSMAP_AD92_Xiong2023 --config_id config02 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset ROSMAP_AD92_Xiong2023 --config_id config03 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset ROSMAP_AD92_Xiong2023 --config_id config04 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset ROSMAP_AD92_Xiong2023 --config_id config05 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset ROSMAP_AD92_Xiong2023 --config_id config06 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset ROSMAP_AD92_Xiong2023 --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset ROSMAP_AD92_Xiong2023 --config_id config07 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset ROSMAP_AD92_Xiong2023 --config_id config07 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset ROSMAP_AD92_Xiong2023 --config_id config07 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset ROSMAP_AD92_Xiong2023 --config_id config07 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset ROSMAP_AD92_Xiong2023 --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset ROSMAP_AD92_Xiong2023 --config_id config07 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset ROSMAP_AD92_Xiong2023 --config_id config07 --n_core 15

Rscript ../DALE_Eval/run/run_InstaPrism.R --dataset ROSMAP_AD92_Xiong2023 --config_id config08 --n_core 15
Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset ROSMAP_AD92_Xiong2023 --config_id config08 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset ROSMAP_AD92_Xiong2023 --config_id config08 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset ROSMAP_AD92_Xiong2023 --config_id config08 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset ROSMAP_AD92_Xiong2023 --config_id config08 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset ROSMAP_AD92_Xiong2023 --config_id config08 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset ROSMAP_AD92_Xiong2023 --config_id config08 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset ROSMAP_AD92_Xiong2023 --config_id config08 --n_core 15

Rscript ../DALE_Eval/run/run_CIBERSORTx.R --dataset ROSMAP_AD92_Xiong2023 --config_id config09 --n_core 15 --extra_args "username=${CIBERSORTX_USERNAME};token=${CIBERSORTX_TOKEN};singularity_container_path=${CIBERSORTX_CONTAINER_PATH}"
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset ROSMAP_AD92_Xiong2023 --config_id config09 --n_core 15 --extra_args 'ENIGMAmode=L2'
Rscript ../DALE_Eval/run/run_ENIGMA.R --dataset ROSMAP_AD92_Xiong2023 --config_id config09 --n_core 15 --extra_args 'ENIGMAmode=trace'
Rscript ../DALE_Eval/run/run_bMIND.R --dataset ROSMAP_AD92_Xiong2023 --config_id config09 --n_core 15 --extra_args 'export_posterior=true;run_epicunmix=false'
Rscript ../DALE_Eval/run/run_TCA.R --dataset ROSMAP_AD92_Xiong2023 --config_id config09 --n_core 15
Rscript ../DALE_Eval/run/run_Unico.R --dataset ROSMAP_AD92_Xiong2023 --config_id config09 --n_core 15
Rscript ../DALE_Eval/run/run_EPICunmix_with_posterior.R --dataset ROSMAP_AD92_Xiong2023 --config_id config09 --n_core 15
