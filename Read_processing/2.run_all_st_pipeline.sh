#!/usr/bin/bash

conda activate stpipeline_1-7-9_visium
sbatch st_pipeline_Visium.sh EXP13_121_A1
sbatch st_pipeline_Visium.sh EXP13_121_A2
sbatch st_pipeline_Visium.sh EXP13_121_B1
sbatch st_pipeline_Visium.sh EXP13_121_B2
sbatch st_pipeline_Visium.sh EXP13_121_C1
sbatch st_pipeline_Visium.sh EXP13_121_C2
sbatch st_pipeline_Visium.sh EXP13_121_D1
sbatch st_pipeline_Visium.sh EXP12_128_A1
sbatch st_pipeline_Visium.sh EXP12_128_A2
sbatch st_pipeline_Visium.sh EXP12_128_B1
sbatch st_pipeline_Visium.sh EXP12_128_B2
sbatch st_pipeline_Visium.sh EXP12_128_C1
sbatch st_pipeline_Visium.sh EXP12_128_C2
sbatch st_pipeline_Visium.sh EXP9_074_A2
sbatch st_pipeline_Visium.sh EXP10_076_B1
sbatch st_pipeline_Visium.sh EXP10_076_C2

