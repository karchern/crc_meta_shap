for r_id in 1 2 3 4 5 
do	
	for f_id in 1 2 3 4 5 6 7 8 9 10
	do
		for on_what in training testing
		do 
			echo -e "Rscript run_SHAP_test.r ${r_id} ${f_id} ${on_what}"
		done
	done
done > iv.sh
