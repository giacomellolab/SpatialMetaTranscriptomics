#!/usr/bin/bash

conda activate stpipeline_1-7-9_visium
sbatch probe_table.sh EXP13_121_A1
sbatch probe_table.sh EXP13_121_A2
sbatch probe_table.sh EXP13_121_B1

sbatch probe_table.sh EXP13_121_B2_1
sbatch probe_table.sh EXP13_121_B2_21
sbatch probe_table.sh EXP13_121_B2_22

sbatch probe_table.sh EXP13_121_C1_1
sbatch probe_table.sh EXP13_121_C1_21
sbatch probe_table.sh EXP13_121_C1_22

sbatch probe_table.sh EXP13_121_C2

sbatch probe_table.sh EXP13_121_D1_1
sbatch probe_table.sh EXP13_121_D1_21
sbatch probe_table.sh EXP13_121_D1_22

sbatch probe_table.sh EXP12_128_A1_1
sbatch probe_table.sh EXP12_128_A1_21
sbatch probe_table.sh EXP12_128_A1_22

sbatch probe_table.sh EXP12_128_A2
sbatch probe_table.sh EXP12_128_B1
sbatch probe_table.sh EXP12_128_B2

sbatch probe_table.sh EXP12_128_C1_1
sbatch probe_table.sh EXP12_128_C1_21
sbatch probe_table.sh EXP12_128_C1_22

sbatch probe_table.sh EXP12_128_C2_1
sbatch probe_table.sh EXP12_128_C2_21
sbatch probe_table.sh EXP12_128_C2_22

sbatch probe_table.sh EXP9_074_A2_1
sbatch probe_table.sh EXP9_074_A2_2

sbatch probe_table.sh EXP10_076_B1_1
sbatch probe_table.sh EXP10_076_B1_2

sbatch probe_table.sh EXP10_076_C2_1
sbatch probe_table.sh EXP10_076_C2_2

