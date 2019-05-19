/home/chen/nucleosome_model/nucleosome_prediction.pl  -t "occupancy_probs" -s input.fa -tab -c $2 -temp $4 -tab -p occ_probs
mv /home/chen/nucleosome_model/occ_probs.tab /home/chen/Desktop/Work_Project
/home/chen/nucleosome_model/nucleosome_prediction.pl  -t "scores" -raw_binding  -s input.fa -tab -c $2 -temp $4 -tab -p scores
mv /home/chen/nucleosome_model/scores.tab /home/chen/Desktop/Work_Project
