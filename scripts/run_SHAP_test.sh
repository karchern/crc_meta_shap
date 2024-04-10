for r_id in 1 2 3 4 5 
do	
	for f_id in 1 2 3 4 5
	do
		for on_what in training testing
			do
			for which_model in RF lasso
				do
				for dataset in Feng Zeller
				do
					echo -e "Rscript run_SHAP_test.r ${r_id} ${f_id} ${on_what} ${which_model} ${dataset}"
				done
			done
		done
	done
done > iv.sh
