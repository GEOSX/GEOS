#!/bin/bash

python3 uberenv.py \
  -j 36 \
  --spack-env-file "spack_envs/lassen-toss_4_ppc64le_ib.yaml" \
  --spec "@develop %clang@14.0.6 ~openmp +cuda cuda_arch=70 \
          ^essl \
          ^chai@develop \
          ^raja@develop \
          ^umpire@develop \
          ^hypre@develop+cuda cuda_arch=70"

# python3 uberenv.py \
#   -j 36 \
#   --spack-env-file "spack_envs/quartz-toss_4_x86_64_ib.yaml" \
#   --spec "@develop %clang@14.0.6 +openmp \
#           ^intel-oneapi-mkl \
#           ^chai@develop \
#           ^raja@develop \
#           ^umpire@develop \
#           ^hypre@develop+openmp"
