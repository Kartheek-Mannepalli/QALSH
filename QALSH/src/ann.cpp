#include "headers.h"


// -----------------------------------------------------------------------------
int ground_truth(					// output the ground truth results
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	char* data_set,						// address of data set
	char* query_set,					// address of query set
	char* truth_set)					// address of ground truth file
{
//	clock_t startTime = (clock_t) -1;
//	clock_t endTime   = (clock_t) -1;

    auto start_ground_truth_time = std::chrono::steady_clock::now(); // Ground truth time (omid)

	int i, j;
	FILE* fp = nullptr;

	// -------------------------------------------------------------------------
	//  Read data set and query set
	// -------------------------------------------------------------------------
//	startTime = clock();
	g_memory += SIZEFLOAT * (n + qn) * d;
	if (check_mem()) return 1;

	float** data = new float*[n];
	for (i = 0; i < n; i++) data[i] = new float[d];
	if (read_set(n, d, data_set, data)) {
	    printf("Reading Dataset Error!\n");
		exit(EXIT_FAILURE);
	}

	float** query = new float*[qn];
	for (i = 0; i < qn; i++) query[i] = new float[d];
	if (read_set(qn, d, query_set, query) == 1) {
	    printf("Reading Query Set Error!\n");
		exit(EXIT_FAILURE);
	}
//	endTime = clock();
//	printf("Read Dataset and Query Set: %.6f Seconds\n\n",
//		((float) endTime - startTime) / CLOCKS_PER_SEC);

	// -------------------------------------------------------------------------
	//  output ground truth results (using linear scan method)
	// -------------------------------------------------------------------------
	int maxk = MAXK;

    if (isQinDS) // whether query is in dataset or not (omid)
        maxk++;

	float dist = -1.0F;
	float* knndist = new float[maxk];
	g_memory += SIZEFLOAT * maxk;

	fp = fopen(truth_set, "w");		// open output file
	if (!fp) {
		printf("I could not create %s.\n", truth_set);
		return 1;
	}

	fprintf(fp, "%d %d\n", qn, maxk);
	for (i = 0; i < qn; i++) {
		for (j = 0; j < maxk; j++) {
			knndist[j] = MAXREAL;
		}
									// find k-nn points of query
		for (j = 0; j < n; j++) {
			dist = calc_l2_dist(data[j], query[i], d);

			int ii, jj;
			for (jj = 0; jj < maxk; jj++) {
				if (compfloats(dist, knndist[jj]) == -1) {
					break;
				}
			}
			if (jj < maxk) {
				for (ii = maxk - 1; ii >= jj + 1; ii--) {
					knndist[ii] = knndist[ii - 1];
				}
				knndist[jj] = dist;
			}
		}

		fprintf(fp, "%d", i + 1);	// output Lp dist of k-nn points
		for (j = 0; j < maxk; j++) {
			fprintf(fp, " %f", knndist[j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);						// close output file
//	endTime = clock();
//	printf("Generate Ground Truth: %.6f Seconds\n\n",
//		((float) endTime - startTime) / CLOCKS_PER_SEC);
    auto end_ground_truth_time = std::chrono::steady_clock::now(); // Ground truth time (omid)
    GT_TIME += std::chrono::duration <double, std::milli>
            (end_ground_truth_time - start_ground_truth_time).count(); // Ground truth time (omid)

    printf("*** Ground truth time:\t%0.6lf\n", GT_TIME); // Ground truth time (omid)


	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	if (data != nullptr) {				// release <data>
		for (i = 0; i < n; i++) {
			delete[] data[i]; data[i] = nullptr;
		}
		delete[] data; data = nullptr;
		g_memory -= SIZEFLOAT * n * d;
	}
	if (query != nullptr) {			// release <query>
		for (i = 0; i < qn; i++) {
			delete[] query[i]; query[i] = nullptr;
		}
		delete[] query; query = nullptr;
		g_memory -= SIZEFLOAT * qn * d;
	}
	if (knndist != nullptr) {			// release <knndist>
		delete[] knndist; knndist = nullptr;
		g_memory -= SIZEFLOAT * maxk;
	}

	//printf("memory = %.2f MB\n", (float) g_memory / (1024.0f * 1024.0f));
	return 0;
}

// -----------------------------------------------------------------------------
int indexing(						// build hash tables for the dataset
	int   n,							// number of data points
	int   d,							// dimension of space
	int   B,							// page size
	float ratio,						// approximation ratio
	char* data_set,						// address of data set
	char* output_folder)				// folder to store info of qalsh
{
//	clock_t startTime = (clock_t) -1;
//	clock_t endTime   = (clock_t) -1;

    auto start_indexing_time = std::chrono::steady_clock::now(); // Indexing time (omid)

    // Initialize timers (omid)
    READ_DS_TIME = 0;
    WRITE_DS_BIN_TIME = 0;
    INIT_PARAMS_TIME = 0;
    INIT_HASH_TIME = 0;
    PROJ_POINTS_TIME = 0;
    SORT_TIME = 0;
    BUILD_WRITE_TREE_TIME = 0;

	// -------------------------------------------------------------------------
	//  Read data set
	// -------------------------------------------------------------------------
//	startTime = clock();
	g_memory += SIZEFLOAT * n * d;	
	if (check_mem()) return 1;

    auto read_ds_time_start = std::chrono::steady_clock::now();
    float** data = new float*[n];
	for (int i = 0; i < n; i++) data[i] = new float[d];
	if (read_set(n, d, data_set, data) == 1) {
	    printf("Reading Dataset Error!\n");
		exit(EXIT_FAILURE);
	}
    auto read_ds_time_end = std::chrono::steady_clock::now();
    READ_DS_TIME += std::chrono::duration <double, std::milli>
            (read_ds_time_end - read_ds_time_start).count();

//	endTime = clock();
//	printf("Read Dataset: %.6f Seconds\n\n",
//		((float) endTime - startTime) / CLOCKS_PER_SEC);

	// -------------------------------------------------------------------------
	//  Write the data set in new format to disk
	// -------------------------------------------------------------------------
//	startTime = clock();
//	write_data_new_form(n, d, B, data, output_folder);
//	endTime = clock();
//	printf("Write Dataset in New Format: %.6f Seconds\n\n",
//		((float) endTime - startTime) / CLOCKS_PER_SEC);

    // Write data set to one binary file (omid) (begin)
    auto t_write_ds_bin_start = std::chrono::steady_clock::now();

    std::ofstream ds_bin_out;
    ds_bin_out.open(DATA_BIN_PATH, std::ios::binary | std::ios::out);

    for (int pid = 0; pid < n; pid++)
        for (int dim = 0; dim < d; dim++)
            ds_bin_out.write((char*) (&data[pid][dim]), SIZEFLOAT);

    auto t_write_ds_bin_end = std::chrono::steady_clock::now();
    WRITE_DS_BIN_TIME += std::chrono::duration <double, std::milli>
            (t_write_ds_bin_end - t_write_ds_bin_start).count();
    // Write data set to one binary file (omid) (end)

	// -------------------------------------------------------------------------
	//  Bulkloading
	// -------------------------------------------------------------------------
//	char fname[200];
//	strcpy(fname, output_folder);
//	strcat(fname, "L2_index.out");
//
//	FILE* fp = fopen(fname, "w");
//	if (!fp) {
//		printf("I could not create %s.\n", fname);
//		return 1;					// fail to return
//	}

//	startTime = clock();
	QALSH* lsh = new QALSH();
	lsh->init(n, d, B, ratio, output_folder);
	lsh->bulkload(data);
//	endTime = clock();

    auto end_indexing_time = std::chrono::steady_clock::now(); // Indexing time (omid)
    INDEXING_TIME += std::chrono::duration <double, std::milli>
            (end_indexing_time - start_indexing_time).count(); // indexing time (omid)

    printf("*** Indexing time:\t%0.6lf\n", INDEXING_TIME); // Indexing time (omid)

    printf("READ_DS_TIME: %.6f\n", READ_DS_TIME);
    printf("WRITE_DS_BIN_TIME: %.6f\n", WRITE_DS_BIN_TIME);
    printf("INIT_PARAMS_TIME: %.6f\n", INIT_PARAMS_TIME);
    printf("INIT_HASH_TIME: %.6f\n", INIT_HASH_TIME);
    printf("PROJ_POINTS_TIME: %.6f\n", PROJ_POINTS_TIME);
    printf("SORT_TIME: %.6f\n", SORT_TIME);
    printf("BUILD_WRITE_TREE_TIME: %.6f\n", BUILD_WRITE_TREE_TIME);

//	float indexing_time = ((float) endTime - startTime) / CLOCKS_PER_SEC;
//	printf("\nIndexing Time: %.6f seconds\n\n", indexing_time);
//	fprintf(fp, "Indexing Time: %.6f seconds\n", indexing_time);
//	fclose(fp);

	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	if (data != nullptr) {
		for (int i = 0; i < n; i++) {
			delete[] data[i]; data[i] = nullptr;
		}
		delete[] data; data = nullptr;
		g_memory -= SIZEFLOAT * n * d;
	}
	if (lsh != nullptr) {
		delete lsh; lsh = nullptr;
	}

	//printf("memory = %.2f MB\n", (float) g_memory / (1024.0f * 1024.0f));
	return 0;
}

// -----------------------------------------------------------------------------
int lshknn(							// k-nn via qalsh (data in disk)
	int   qn,							// number of query points
	int   d,							// dimensionality
	char* query_set,					// path of query set
	char* truth_set,					// groundtrue file
	char* output_folder)				// output folder
{
	int ret = 0;
	int maxk = MAXK;

    if (isQinDS) // whether query is in dataset or not (omid)
        maxk++;

	int i, j;
	FILE* fp = nullptr;				// file pointer

	// -------------------------------------------------------------------------
	//  Read query set
	// -------------------------------------------------------------------------
	g_memory += SIZEFLOAT * qn * d;
	float** query = new float*[qn];
	for (i = 0; i < qn; i++) query[i] = new float[d];
	if (read_set(qn, d, query_set, query)) {
	    printf("Reading Query Set Error!\n");
		exit(EXIT_FAILURE);
	}

	// -------------------------------------------------------------------------
	//  Read the ground truth file
	// -------------------------------------------------------------------------
	g_memory += SIZEFLOAT * qn * maxk;
//	float* R = new float[qn * maxk];

	fp = fopen(truth_set, "r");		// open ground truth file
	if (!fp) {
		printf("Could not open the ground truth file.\n");
		return 1;
	}

//	fscanf(fp, "%d %d\n", &qn, &maxk);
//	for (int i = 0; i < qn; i++) {
//		fscanf(fp, "%d", &j);
//		for (j = 0; j < maxk; j ++) {
//			fscanf(fp, " %f", &(R[i * maxk + j]));
//		}
//	}

    int temp_qn; // omid
    fscanf(fp, "%d %d\n", &temp_qn, &maxk);
    float* R = new float[temp_qn * maxk];
    for (int i = 0; i < temp_qn; i++) {
        fscanf(fp, "%d", &j);
        for (j = 0; j < maxk; j ++) {
            fscanf(fp, " %f", &(R[i * maxk + j]));
        }
    }

	fclose(fp);						// close groundtrue file

	// -------------------------------------------------------------------------
	//  K-nearest neighbor (k-nn) search via qalsh
	// -------------------------------------------------------------------------
	int kNNs[] = {1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
	int maxRound = 11;
	int top_k = 0;

	float allTime   = -1.0f;
	float allRatio  = -1.0f;
	float thisRatio = -1.0f;
	float knnStartTime = -1.0f;
	float knnEndTime = -1.0f;
	float knnTime = 0.0f;

    FILE* exact_fp = nullptr; // omid
    int exact_results[101]; // id (omid)
    char gt_dir[100] = ""; // omid
    get_leading_folder(truth_set, gt_dir); // omid
    float map = -1; // omid

//	clock_t startTime = (clock_t) -1.0;
//	clock_t endTime   = (clock_t) -1.0;
									// init the results
	g_memory += (long) sizeof(ResultItem) * maxk;
	ResultItem* rslt = new ResultItem[maxk];
	for (i = 0; i < maxk; i++) {
		rslt[i].id_ = -1;
		rslt[i].dist_ = MAXREAL;
	}

	QALSH* lsh = new QALSH();		// restore QALSH
	if (lsh->restore(output_folder)) {
	    printf("Could not restore qalsh\n");
		exit(EXIT_FAILURE);
	}

//	char output_set[200];
//	strcpy(output_set, output_folder);
//	strcat(output_set, "L2_qalsh.out");
//
//	fp = fopen(output_set, "w");	// open output file
//	if (!fp) {
//		printf("Could not create the output file.\n");
//		return 1;
//	}

    // Use QALSH output format (omid) (begin)
    printf("TOP_K\tRATIO\tmAP\tQUERY_TIME(ms)\tALG_TIME(ms)\t"
           "INDEX_IO_TIME(ms)\tINDEX_IO_SIZE(KB)\tINDEX_IO_NUM\t"
           "DATA_IO_TIME(ms)\tDATA_IO_SIZE(KB)\tDATA_IO_NUM\t"
           "DIST_CALC_TIME\t"
           "VR_ALG_TIME\tVR_IO_TIME\tVR_TIME\tVR_COUNT\n");
    // Use QALSH output format (omid) (end)

    // Excel-ready output (omid) (begin)
    char excel_output[200];
    strcpy(excel_output, output_folder);

    if (excel_output[strlen(excel_output) - 1] == '/')
        strcat(excel_output, "excel.out");
    else
        strcat(excel_output, "/excel.out");

    FILE* excel_file = fopen(excel_output, "w");

    if (!excel_file) {
        std::cout << "Can't open " << excel_output << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string ds_name(query_set);
    ds_name = ds_name.substr(ds_name.rfind("/") + 1, ds_name.rfind(".") - ds_name.rfind("/") - 1);
    fprintf(excel_file, "QALSH_%s\n", ds_name.c_str());

    fprintf(excel_file,"TOP_K\tRATIO\tmAP\tQUERY_TIME(ms)\tALG_TIME(ms)\t"
                       "INDEX_IO_TIME(ms)\tINDEX_IO_SIZE(KB)\tINDEX_IO_NUM\t"
                       "DATA_IO_TIME(ms)\tDATA_IO_SIZE(KB)\tDATA_IO_NUM\t"
                       "DIST_CALC_TIME\t"
                       "VR_ALG_TIME\tVR_IO_TIME\tVR_TIME\tVR_COUNT\n");
    // Excel-ready output (omid) (end)


//	printf("QALSH for c-k-ANN Search: \n");
//	printf("    Top-k\tRatio\t\tI/O\tDistIO\tPageIO\t\tTime (ms)\tReadDist (ms)\tCalcL2 (ms)\tKNN (ms)\n");

	for (int num = 0; num < maxRound; num++) {
		auto output_k = top_k = kNNs[num];

        // Fix ratio calculation if query is inside the dataset (omid)
        if (isQinDS)
            top_k++;

		allRatio = 0.0f;
		map = 0.0f; // omid
		knnTime = 0;
//		startTime = clock();

        // Initialize timers (begin) (omid)
        QUERY_TIME = 0;
        ALG_TIME = 0;
        DIST_CALC_TIME = 0;
        INDEX_IO_TIME = 0;
        DATA_IO_TIME = 0;
        INDEX_IO_NUM = 0;
        DATA_IO_NUM = 0;
        INDEX_IO_SIZE = 0;
        DATA_IO_SIZE = 0;
        VR_ALG_TIME = 0;
        VR_IO_TIME = 0;
        VR_TIME = 0;
        VR_COUNT = 0;
        VR_FLAG = false;
        // Initialize timers (end) (omid)

		for (i = 0; i < qn; i++) {
//			knnStartTime = clock();

            auto start_query_time = std::chrono::steady_clock::now(); // Total query time (omid)

			lsh->knn(query[i], top_k, rslt, output_folder);

            auto end_query_time = std::chrono::steady_clock::now(); // Total query time (omid)
            QUERY_TIME += std::chrono::duration <double, std::milli>
                    (end_query_time - start_query_time).count(); // Total query time (omid)

//			knnEndTime = clock();
//			knnTime += (((float) knnEndTime - knnStartTime) / CLOCKS_PER_SEC) * 1000.0f;
			thisRatio = 0.0f;

            int zeroCount = 0; // omid
            for (j = 0; j < top_k; j++) {
                // omid
                if (R[i * maxk + j] == 0)
                    zeroCount++;
                else
                    // omid
                    thisRatio += rslt[j].dist_ / R[i * maxk + j];
            }
//			thisRatio /= top_k;
            thisRatio /= (top_k - zeroCount); // Keep query in the ds (omid)
			allRatio += thisRatio;
			//printf("\nALL IO = %d\n",allIO);
			//printf("        No.%2d: Top-k = %d, IO = %4d, Ratio = %0.6f\n",
			//	i+1, top_k, thisIO, thisRatio);

            // Calc mAP (omid) (begin)
            char temp_path[100];
            char temp_buff[50];
            strcpy(temp_path, gt_dir);
            sprintf(temp_buff, "exact_results/%d", qn);
            strcat(temp_path, temp_buff);

            exact_fp = fopen(temp_path, "r");
            if (!exact_fp) {
                printf("Could not open the file of exact results.\n");
                exit(EXIT_FAILURE);
            }

            for (i = 0; i < maxk; i++) {
                if (feof(exact_fp)) {
                    printf("Exact results file is incomplete\n");
                    exit(EXIT_FAILURE);
                }

                fscanf(exact_fp, "%d", &exact_results[i]);
            }

            float ap_score = 0.0f;
            float ap_num_hits = 0.0f;

            for (i = 0; i < top_k; i++) {
                for (j = 0; j < top_k; j++) {
                    if (rslt[i].id_ == exact_results[j]) {
                        ap_num_hits += 1.0f;
                        ap_score += ap_num_hits / (i + 1.0f);
                    }
                }
            }

            ap_score /= top_k;
            map += ap_score;

            fclose(exact_fp);
            // Calc mAP (omid) (end)
		}

//		endTime = clock();
//		allTime = ((float) endTime - startTime) / CLOCKS_PER_SEC;

		//average the times to a per query basis
		allRatio = allRatio / qn;
        map /= qn; // omid
//		allTime = (allTime * 1000.0f) / qn;
//		knnTime = knnTime/qn;
//		printf("    %3d\t\t%.4f\t\t%d\t%d\t%d\t\t%.2f\t\t%lf\t%lf\t%lf\n", top_k, allRatio,
//			allIO, dist_io, page_io, allTime, read_dataTime, calc_l2_distTime,knnTime);
//		fprintf(fp, "%3d\t\t%.4f\t\t%d\t%d\t%d\t\t%.2f\t\t%lf\t%lf\t%lf\n", top_k, allRatio,
//				allIO, dist_io, page_io, allTime, read_dataTime, calc_l2_distTime,knnTime);

        VR_TIME = VR_ALG_TIME + VR_IO_TIME; // omid

        printf("%d\t%0.6lf\t%0.6lf\t%0.6lf\t%0.6lf\t"
               "%0.6lf\t%0.2lf\t%d\t"
               "%0.6lf\t%0.2lf\t%d\t"
               "%0.6lf\t"
               "%0.6lf\t%0.6lf\t%0.6lf\t%d\n",
               output_k, allRatio, map, QUERY_TIME, ALG_TIME,
               INDEX_IO_TIME, INDEX_IO_SIZE, INDEX_IO_NUM,
               DATA_IO_TIME, DATA_IO_SIZE, DATA_IO_NUM,
               DIST_CALC_TIME,
               VR_ALG_TIME, VR_IO_TIME, VR_TIME, VR_COUNT);

        // Excel-ready output (omid) (begin)
        fprintf(excel_file,"%d\t%0.6lf\t%0.6lf\t%0.6lf\t%0.6lf\t"
                           "%0.6lf\t%0.2lf\t%d\t"
                           "%0.6lf\t%0.2lf\t%d\t"
                           "%0.6lf\t"
                           "%0.6lf\t%0.6lf\t%0.6lf\t%d\n",
                output_k, allRatio, map, QUERY_TIME, ALG_TIME,
                INDEX_IO_TIME, INDEX_IO_SIZE, INDEX_IO_NUM,
                DATA_IO_TIME, DATA_IO_SIZE, DATA_IO_NUM,
                DIST_CALC_TIME,
                VR_ALG_TIME, VR_IO_TIME, VR_TIME, VR_COUNT);
        // Excel-ready output (omid) (end)

	}

	printf("\n");
//	fclose(fp);						// close output file

    fclose(excel_file); // omid
    if (lsh->write_para_to_excel(excel_output)) return 1; // omid

	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	if (query != nullptr) {			// release <query>
		for (i = 0; i < qn; i++) {
			delete[] query[i]; query[i] = nullptr;
		}
		delete[] query; query = nullptr;
		g_memory -= SIZEFLOAT * qn * d;
	}
	if (lsh != nullptr) {				// release <lsh>
		delete lsh; lsh = nullptr;
	}
									// release <R> and (/or) <rslt>
	if (R != nullptr || rslt != nullptr) {
		delete[] R; R = nullptr;
		delete[] rslt; rslt = nullptr;
		g_memory -= (SIZEFLOAT * qn * maxk + sizeof(ResultItem) * maxk);
	}

	//printf("memory = %.2f MB\n", (float) g_memory / (1024.0f * 1024.0f));
	return ret;
}

// -----------------------------------------------------------------------------
int linear_scan(					// brute-force linear scan (data in disk)
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   B,							// page size
	char* query_set,					// address of query set
	char* truth_set,					// address of ground truth file
    char* data_set,						// address of data set
	char* output_folder)				// output folder
{
	// -------------------------------------------------------------------------
	//  Allocation and initialzation.
	// -------------------------------------------------------------------------
	clock_t startTime = (clock_t) -1.0f;
	clock_t endTime   = (clock_t) -1.0f;

	int kNNs[] = {1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
	int maxRound = 11;

	int i, j, top_k;
	int maxk = MAXK;

    if (isQinDS) // whether query is in dataset or not (omid)
        maxk++;

	float allTime   = -1.0f;
	float thisRatio = -1.0f;
	float allRatio  = -1.0f;

	g_memory += (SIZEFLOAT * (d + d + (qn + 1) * maxk) + SIZECHAR * (600 + B));
	
	float* knn_dist = new float[maxk];
	for (i = 0; i < maxk; i++) {
		knn_dist[i] = MAXREAL;
	}

	float** R = new float*[qn];
	for (i = 0; i < qn; i++) {
		R[i] = new float[maxk];
		for (j = 0; j < maxk; j++) {
			R[i][j] = 0.0f;
		}
	}

	float* data     = new float[d];	// one data object
	float* query    = new float[d];	// one query object

	char* buffer    = new char[B];	// every time can read one page
	char* fname     = new char[200];// file name for data
	char* data_path = new char[200];// data path
	char* out_set	= new char[200];// output file

    auto** temp_data = new float*[n];
    for (int temp_i = 0; temp_i < n; temp_i++) temp_data[temp_i] = new float[d];
    if (read_set(n, d, data_set, temp_data) == 1) {
        printf("Reading Dataset Error!\n");
        exit(EXIT_FAILURE);
    }
    write_data_new_form(n, d, B, temp_data, output_folder);

    // -------------------------------------------------------------------------
	//  Open the output file, and read the ground true results
	// -------------------------------------------------------------------------
//	strcpy(out_set, output_folder);	// generate output file
//	strcat(out_set, "L2_linear.out");
//
//	FILE* ofp = fopen(out_set, "w");
//	if (!ofp) {
//		printf("I could not create %s.\n", out_set);
//		return 1;
//	}
									// open ground true file
	FILE* tfp = fopen(truth_set, "r");
	if (!tfp) {
		printf("I could not create %s.\n", truth_set);
		return 1;
	}
									// read top-k nearest distance
	fscanf(tfp, "%d %d\n", &qn, &maxk);
	for (int i = 0; i < qn; i++) {
		fscanf(tfp, "%d", &j);
		for (j = 0; j < maxk; j ++) {
			fscanf(tfp, " %f", &(R[i][j]));
		}
	}
	fclose(tfp);					// close ground true file

	// -------------------------------------------------------------------------
	//  Calc the number of data object in one page and the number of data file.
	//  <num> is the number of data in one data file
	//  <total_file> is the total number of data file
	// -------------------------------------------------------------------------
	int num = (int) floor((float) B / (d * SIZEFLOAT));
	int total_file = (int) ceil((float) n / num);
	if (total_file == 0) return 1;

	// -------------------------------------------------------------------------
	//  Brute-force linear scan method (data in disk)
	//  For each query, we limit that we can ONLY read one page of data.
	// -------------------------------------------------------------------------
	int count = 0;
	float dist = -1.0F;
									// generate the data path
	strcpy(data_path, output_folder);
	strcat(data_path, "data/");

	printf("Linear Scan Search:\n");
	printf("    Top-k\tRatio\t\tI/O\t\tTime (ms)\n");
	for (int round = 0; round < maxRound; round++) {
		top_k = kNNs[round];
		allRatio = 0.0f;

//		startTime = clock();
		FILE* qfp = fopen(query_set, "r");
		if (!qfp) {
		    printf("Could not open the query set.\n");
		    exit(EXIT_FAILURE);
		}

		for (i = 0; i < qn; i++) {
			// -----------------------------------------------------------------
			//  Step 1: read a query from disk and init the k-nn results
			// -----------------------------------------------------------------
			fscanf(qfp, "%d", &j);
			for (j = 0; j < d; j++) {
				fscanf(qfp, " %f", &query[j]);
			}

			for (j = 0; j < top_k; j++) {
				knn_dist[j] = MAXREAL;
			}

			// -----------------------------------------------------------------
			//  Step 2: find k-nn results for the query
			// -----------------------------------------------------------------
			for (j = 0; j < total_file; j++) {
				// -------------------------------------------------------------
				//  Step 2.1: get the file name of current data page
				// -------------------------------------------------------------
				get_data_filename(j, data_path, fname);

				// -------------------------------------------------------------
				//  Step 2.2: read one page of data into buffer
				// -------------------------------------------------------------
				if (read_buffer_from_page(B, fname, buffer) == 1) {
				    printf("error to read a data page\n");
					exit(EXIT_FAILURE);
				}

				// -------------------------------------------------------------
				//  Step 2.3: find the k-nn results in this page. NOTE: the 
				// 	number of data in the last page may be less than <num>
				// -------------------------------------------------------------
				if (j < total_file - 1) count = num;
				else count = n % num;

				for (int z = 0; z < count; z++) {
					read_data_from_buffer(z, d, data, buffer);
					dist = calc_l2_dist(data, query, d);

					int ii, jj;
					for (jj = 0; jj < top_k; jj++) {
						if (compfloats(dist, knn_dist[jj]) == -1) {
							break;
						}
					}
					if (jj < top_k) {
						for (ii = top_k - 1; ii >= jj + 1; ii--) {
							knn_dist[ii] = knn_dist[ii - 1];
						}
						knn_dist[jj] = dist;
					}
				}
			}

			thisRatio = 0.0f;

            int zeroCount = 0; // Keep query in the ds (omid)
            for (j = 0; j < top_k; j++) {
                // Keep query in the ds (omid) (begin)
                if (R[i * maxk + j] == 0)
                    zeroCount++;
                else
                    // Keep query in the ds (omid) (end)
                    thisRatio += knn_dist[j] / R[i][j];
            }
//			thisRatio /= top_k;
            thisRatio /= (top_k - zeroCount); // Keep query in the ds (omid)

			allRatio += thisRatio;
		}
		// -----------------------------------------------------------------
		//  Step 3: output result of top-k nn points
		// -----------------------------------------------------------------
		fclose(qfp);				// close query file
//		endTime  = clock();
		allTime  = ((float) endTime - startTime) / CLOCKS_PER_SEC;
		allTime = (allTime * 1000.0f) / qn;
		allRatio = allRatio / qn;
									// output results
		printf("    %3d\t\t%.4f\t\t%d\t\t%.2f\n", top_k, allRatio, 
			total_file, allTime);
//		fprintf(ofp, "%d\t%f\t%d\t%f\n", top_k, allRatio, total_file, allTime);
	}
	printf("\n");
//	fclose(ofp);					// close output file

	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	if (R != nullptr) {
		for (i = 0; i < qn; i++) {
			delete[] R[i]; R[i] = nullptr;
		}
		delete[] R; R = nullptr;
	}
	if (knn_dist != nullptr || buffer != nullptr || data != nullptr || query != nullptr) {
		delete[] knn_dist; knn_dist = nullptr;
		delete[] buffer; buffer = nullptr;
		delete[] data; data = nullptr;
		delete[] query; query = nullptr;
	}
	if (fname != nullptr || data_path != nullptr || out_set != nullptr) {
		delete[] fname; fname = nullptr;
		delete[] data_path; data_path = nullptr;
		delete[] out_set; out_set = nullptr;
	}
	g_memory -= (SIZEFLOAT * (d + d + (qn + 1) * maxk) + SIZECHAR * (600 + B));
	
	//printf("memory = %.2f MB\n", (float) g_memory / (1024.0f * 1024.0f));
	return 0;
}
