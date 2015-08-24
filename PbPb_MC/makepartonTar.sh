#!bin/bash

root -l <<EOF
.L RAA_partonread_mc_pbpb.C++
.q
EOF

rm jetRAA_run_PbPbparton_MC.tar
tar -zcvf jetRAA_run_PbPbparton_MC.tar jetRAA_PbPb_*.txt RAA_partonread_*.* weights.root data_mc_cent_vz_weight.root weights_pbpb_*.root
