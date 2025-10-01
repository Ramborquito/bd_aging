/* Elaborado por Marco Antonio Ramirez Guizar
Codigo en C++ paralelizado con nvidia-CUDA para simular un sistema bidisperso de particulas
para estudiar las fuerzas efectivas bajo el esquema de contraccion de fuerzas en un sistema
en enfriamiento o calentamiento.
*/
# include "encabezados.h"

# define LOG_SAMPLE
# define GRO_FLAG 0
# define HOST 0

//===================================== 2025-03-11 ===================================

std::vector<float> log_sampling(int n) {
    std::vector<float> values;
    values.reserve(n + 1);

    float r = std::pow(10.0f, 1.0f / n); // raíz n-ésima de 10
    float x = 1.0f;

    for (int i = 0; i <= n; ++i) {
        values.push_back(x);
        x *= r;
    }

    return values;
}

std::vector<int> sampling_indexes(float minimal_decade, float run_time, float dt, int samplings_per_decade) {
    int last_power = static_cast<int>(round(log10(run_time)));
    int initial_power = static_cast<int>(round(log10(minimal_decade)));
    const int decades = last_power - initial_power;

    auto values_per_scale = log_sampling(samplings_per_decade);

    int total_indexes = decades * samplings_per_decade;

    std::vector<int> sampling_indexes(total_indexes + 1);
    int sampling_counter = 0;
    sampling_indexes[0] = 0; // initial configuration
    for (int ipow = initial_power; ipow < last_power; ++ipow) {
        float power = pow(10, ipow);
        for (int j = 0; j < samplings_per_decade; ++j) {
            float scaled_time;
            int index;
            scaled_time = values_per_scale[j] * power;
            index = (int) round(scaled_time / dt);
            sampling_counter++;
            sampling_indexes[sampling_counter] = index;
        }
    }
    return sampling_indexes;
}


int main() {
    float3 *rr_big_vec, *rr_big_raw_vec, *rr_big_ini_vec, *ff_big_vec,
            *rr_sml_vec, *rr_sml_raw_vec, *rr_sml_ini_vec, *ff_sml_vec;
    float3 rr, drr, zero;
    float *vir_big_vec, *vir_sml_vec, *pot_big_vec, *pot_sml_vec;
    float *gder_bb_vec, *gder_bs_vec, *gder_sb_vec, *gder_ss_vec;
    int *nocup_big_vec, *cell_big_vec, *nocup_sml_vec, *cell_sml_vec;
    double run_time, trans_time, dt;
    float side, cutoff, bin_size_gder, qmax,
            shell_vol, vol_free, msd_big,
            msd_sml, dist, aux, sigma, sigma_big, sigma_sml, mass_big, mass_sml,
            range_gder, xngrain_big, xngrain_sml, xngrain_tot, big_z,
            virial, ene_pot_big, ene_pot_sml, xnb, T0, Tf, t0, tf, time, period,
            cell_side_big, phi_big, phi_sml, cell_side_sml, volume;
    float ***gders, *energy_ub, *energy_us, *energy_temp, *msd_b, *msd_s,
            *time_energy, *time_msd;
    int ngrain_tot, ngrain_big, ngrain_sml, ni, niter, ntrans, ngap, idum,
            nsamples, ncell_big, ncell_sml, ncell_big3, ncell_sml3, ntot_big,
            ntot_sml, ii, jj, nbins_gder, ntags_big,
            ntags_sml, nocup, nocup_big_max, nocup_sml_max, counter, ngap_rescaling, n_configs;
    int NH, NB_BIG, NB_SML;
    double minimal_decade;
    bool should_sample, print_decades_configs;
    int number_of_tws, number_of_tw, current_decade;

# if !HOST
    int NB_CELL_BIG3, NB_CELL_SML3, NB_NTOT_BIG, NB_NTOT_SML;
# endif
    parametros pars;
    char renglon[200], infile[80], gder_fn[80], energy_fn[80], msd_fn[80], snapshots_fn[80],
            press_fn[80], fself_fn[80];
    char temperature_protocol[40], sample_type[40];
    FILE *fp_bitac, *fp_snaps, *fp_energ, *fp_gder, *fp_msd, *fp_in, *fp_out,
            *fp_press, *fp_data, *fp_colors, *fp_fself;
#ifdef LOG_SAMPLE
    fp_data = fopen("wca_aging_log.data", "r");
#else
    fp_data = fopen("wca_aging_linear.data", "r");
#endif
    if (fp_data == nullptr) fp_data = stdin;

    printf("temp_protocol ?\n");
    fgets(renglon, sizeof(renglon), fp_data);
    sscanf(renglon, "%s", temperature_protocol);

    printf("temp_settings ?\n");
    fgets(renglon, sizeof(renglon), fp_data);

    if (strcmp(temperature_protocol, "linear") == 0)
        sscanf(renglon, "%f  %f  %f  %f  %d", &T0, &Tf, &t0, &tf, &ngap_rescaling);

    if (strcmp(temperature_protocol, "sine") == 0)
        sscanf(renglon, "%f  %f  %f  %f  %f  %d", &T0, &Tf, &t0, &tf, &period, &ngap_rescaling);

    printf("dt ?\n");
    fgets(renglon, sizeof(renglon), fp_data);
    sscanf(renglon, "%lf", &dt);

    printf("sample type (linear, log) ?\n");
    fgets(renglon, sizeof(renglon), fp_data);
    sscanf(renglon, "%s", sample_type);

    if (strcmp(sample_type, "log") == 0) {
        printf("log sampling selected. Must be added in compilation as -DLOG_SAMPLE\n");
        printf("trans_time, minimal_decade, run_time, nsamples (taken in run)\n");
        fgets(renglon, sizeof(renglon), fp_data);
        sscanf(renglon, "%lf  %lf %lf", &trans_time, &minimal_decade, &run_time);

        printf("print decade configs?\n");
        fgets(renglon, sizeof(renglon), fp_data);
        sscanf(renglon, "%d  %d  %d", &print_decades_configs, &number_of_tws, &current_decade);
        // is not correct but it works
    } else {
        printf("trans_time, run_time, nsamples (taken in run)\n");
        fgets(renglon, sizeof(renglon), fp_data);
        sscanf(renglon, "%lf %lf %d", &trans_time, &run_time, &nsamples);
    }

    printf("ntags: big, sml ?\n");
    fgets(renglon, sizeof(renglon), fp_data);
    sscanf(renglon, "%d %d", &ntags_big, &ntags_sml);

    printf("n_configs ?\n");
    fgets(renglon, sizeof(renglon), fp_data);
    sscanf(renglon, "%d", &n_configs);

    printf("nbins_gder, range_gder ?\n");
    fgets(renglon, sizeof(renglon), fp_data);
    sscanf(renglon, "%d %f", &nbins_gder, &range_gder);

    printf("qmax ?\n");
    fgets(renglon, sizeof(renglon), fp_data);
    sscanf(renglon, "%f", &qmax);

    printf("idum ?\n");
    fgets(renglon, sizeof(renglon), fp_data);
    sscanf(renglon, "%d", &idum);

    fclose(fp_data);
    printf("Starts simulation_3D_2SP Device\n");

    // calcula y almacena constantes

    pars.dt = dt;
    pars.temp_set = T0;
    pars.nbins_gder = nbins_gder;
    pars.range_gder = range_gder;
    pars.ntags_big = ntags_big;
    pars.ntags_sml = ntags_sml;

    bin_size_gder = range_gder / ((float) nbins_gder);
    pars.bin_size_gder = bin_size_gder;

    zero.x = zero.y = zero.z = 0.0;


    // create a set of idexes for logarithmic sample
#ifdef LOG_SAMPLE
    int points_per_decade = 10;
    auto sampling_idx = sampling_indexes(minimal_decade, run_time, dt, points_per_decade);
    nsamples = (int) sampling_idx.size();
#endif


    //Initialize memory to accumulate results

    //Values to calculate size of arrays
    niter = (int) (run_time / dt);
    ngap = niter / nsamples;

    //g_ij(r) bidisperse
    gders = (float ***) malloc(nsamples * sizeof(float **));
    for (int i = 0; i < nsamples; ++i) {
        gders[i] = (float **) malloc(4 * sizeof(float *));
        for (int j = 0; j < 4; ++j) {
            gders[i][j] = (float *) malloc(nbins_gder * sizeof(float));
            for (int k = 0; k < nbins_gder; ++k) gders[i][j][k] = 0.0;
        }
    }

    //fself
    auto *fself_big = (float *) calloc(nsamples, sizeof(float));
    auto *fself_sml = (float *) calloc(nsamples, sizeof(float));

    //energy and msd and pressure
    time_energy = (float *) malloc(nsamples * sizeof(float));
    energy_ub = (float *) calloc(nsamples, sizeof(float));
    energy_us = (float *) calloc(nsamples, sizeof(float));
    energy_temp = (float *) calloc(nsamples, sizeof(float));

    time_msd = (float *) malloc(nsamples * sizeof(float));
    msd_b = (float *) calloc(nsamples, sizeof(float));
    msd_s = (float *) calloc(nsamples, sizeof(float));

    auto time_pressure = (float *) malloc(nsamples * sizeof(float));
    auto pressure = (float *) calloc(nsamples, sizeof(float));

#if GRO_FLAG
    fp_colors = fopen("color_statistics.out", "w");
#endif
    ////////////////////////////Starts loop//////////////////////////////////////////////////////////////////////

    for (int i_config = 1; i_config <= n_configs; ++i_config) {
        counter = 0;
        int energy_counter = 0;
        number_of_tw = 1;

        // lee datos iniciales
        if (strcmp(sample_type, "log") == 0)
            sprintf(infile, "../configs/decade%d/init_config_%d", current_decade, i_config);
        else
            sprintf(infile, "../configs/init_config_%d", i_config);

        fp_in = fopen(infile, "r");
        if (fp_in == nullptr) {
            printf("Verify file path: %s", infile);
            fflush(stdout);
            exit(-2);
        }

        fgets(renglon, sizeof(renglon), fp_in);
        sscanf(renglon, "%f %f", &phi_big, &phi_sml);
        fgets(renglon, sizeof(renglon), fp_in);
        sscanf(renglon, "%f", &side);
        fgets(renglon, sizeof(renglon), fp_in);
        sscanf(renglon, "%d %d %d", &ngrain_big, &ngrain_sml, &ngrain_tot);
        fgets(renglon, sizeof(renglon), fp_in);
        sscanf(renglon, "%f %f", &sigma_big, &sigma_sml);
        fgets(renglon, sizeof(renglon), fp_in);
        sscanf(renglon, "%f %f", &mass_big, &mass_sml);
        fgets(renglon, sizeof(renglon), fp_in);

        // almacena datos iniciales
        if (i_config == 1) {
            pars.side = side;
            pars.ngrain_big = ngrain_big;
            pars.ngrain_sml = ngrain_sml;
            pars.sigma_big = sigma_big;
            pars.sigma_sml = sigma_sml;
            pars.mass_big = mass_big;
            pars.mass_sml = mass_sml;
            xngrain_big = (float) ngrain_big;
            xngrain_sml = (float) ngrain_sml;
            xngrain_tot = (float) ngrain_tot;
            volume = side * side * side;

            fp_bitac = fopen("bitacora", "w");

            // cell parameters

            cutoff = sigma_big;
            ncell_big = (int) (side / (1.05 * cutoff));
            cell_side_big = side / ((float) ncell_big);
            pars.ncell_big = ncell_big;
            pars.cell_side_big = cell_side_big;

            cutoff = sigma_sml;
            ncell_sml = (int) (side / (1.05 * cutoff));
            cell_side_sml = side / ((float) ncell_sml);
            pars.ncell_sml = ncell_sml;
            pars.cell_side_sml = cell_side_sml;

            fprintf(fp_bitac, "side  %f  ngr b s tot  %d  %d  %d\n", side, ngrain_big,
                    ngrain_sml, ngrain_tot);
            fprintf(fp_bitac, "ncell b, s  %d  %d\n", ncell_big, ncell_sml);
            fprintf(fp_bitac, "cell_side b, s  %f  %f\n", cell_side_big, cell_side_sml);
            fprintf(fp_bitac, "ntags b, s  %d  %d\n", ntags_big, ntags_sml);
            fflush(fp_bitac);

            // ranges

            pars.nrange_bb = 1;
            aux = (sigma_big + sigma_sml) / 2.0;
            pars.nrange_bs = 1 + (int) (aux / cell_side_sml);
            pars.nrange_sb = 1;
            pars.nrange_ss = 1;

            fprintf(fp_bitac, "ranges bb, bs  %d  %d\n", pars.nrange_bb, pars.nrange_bs);
            fprintf(fp_bitac, "ranges sb, ss  %d  %d\n", pars.nrange_sb, pars.nrange_ss);

            // parameters

            ncell_big3 = ncell_big * ncell_big * ncell_big;
            ncell_sml3 = ncell_sml * ncell_sml * ncell_sml;
            ntot_big = ncell_big3 * ntags_big;
            ntot_sml = ncell_sml3 * ntags_sml;
            fprintf(fp_bitac, "ncell_big3  %d  ntot_big  %d\n", ncell_big3, ntot_big);
            fprintf(fp_bitac, "ncell_sml3  %d  ntot_sml  %d\n", ncell_sml3, ntot_sml);

            // memory allocation
#if !HOST

            cudaMallocManaged(&rr_big_vec, ngrain_big * sizeof(float3));
            cudaMallocManaged(&rr_big_raw_vec, ngrain_big * sizeof(float3));
            cudaMallocManaged(&rr_big_ini_vec, ngrain_big * sizeof(float3));
            cudaMallocManaged(&ff_big_vec, ngrain_big * sizeof(float3));
            cudaMallocManaged(&vir_big_vec, ngrain_big * sizeof(float));
            cudaMallocManaged(&pot_big_vec, ngrain_big * sizeof(float));

            cudaMallocManaged(&nocup_big_vec, ncell_big3 * sizeof(int));
            cudaMallocManaged(&cell_big_vec, ntot_big * sizeof(int));

            cudaMallocManaged(&rr_sml_vec, ngrain_sml * sizeof(float3));
            cudaMallocManaged(&rr_sml_raw_vec, ngrain_sml * sizeof(float3));
            cudaMallocManaged(&rr_sml_ini_vec, ngrain_sml * sizeof(float3));
            cudaMallocManaged(&ff_sml_vec, ngrain_sml * sizeof(float3));
            cudaMallocManaged(&vir_sml_vec, ngrain_sml * sizeof(float));
            cudaMallocManaged(&pot_sml_vec, ngrain_sml * sizeof(float));

            cudaMallocManaged(&nocup_sml_vec, ncell_sml3 * sizeof(int));
            cudaMallocManaged(&cell_sml_vec, ntot_sml * sizeof(int));
#else
            force_ls_vec = (float *) malloc(nbins_dplt * sizeof(float));

            rr_big_vec = (float3 *) malloc(ngrain_big * sizeof(float3));
            rr_big_raw_vec = (float3 *) malloc(ngrain_big * sizeof(float3));
            rr_big_ini_vec = (float3 *) malloc(ngrain_big * sizeof(float3));
            ff_big_vec = (float3 *) malloc(ngrain_big * sizeof(float3));
            vir_big_vec = (float *) malloc(ngrain_big * sizeof(float));
            pot_big_vec = (float *) malloc(ngrain_big * sizeof(float));

            nocup_big_vec = (int *) malloc(ncell_big3 * sizeof(int));
            cell_big_vec = (int *) malloc(ntot_big * sizeof(int));

            rr_sml_vec = (float3 *) malloc(ngrain_sml * sizeof(float3));
            rr_sml_raw_vec = (float3 *) malloc(ngrain_sml * sizeof(float3));
            rr_sml_ini_vec = (float3 *) malloc(ngrain_sml * sizeof(float3));
            ff_sml_vec = (float3 *) malloc(ngrain_sml * sizeof(float3));
            vir_sml_vec = (float *) malloc(ngrain_sml * sizeof(float));
            pot_sml_vec = (float *) malloc(ngrain_sml * sizeof(float));

            nocup_sml_vec = (int *) malloc(ncell_sml3 * sizeof(int));
            cell_sml_vec = (int *) malloc(ntot_sml * sizeof(int));
#endif

            // blocks and threads

            NH = 256;
            pars.NH = NH;

            NB_BIG = 1 + (ngrain_big - 1) / NH;
            NB_SML = 1 + (ngrain_sml - 1) / NH;
#if !HOST
            NB_CELL_BIG3 = 1 + (ncell_big3 - 1) / NH;
            NB_CELL_SML3 = 1 + (ncell_sml3 - 1) / NH;
            NB_NTOT_BIG = 1 + (ntot_big - 1) / NH;
            NB_NTOT_SML = 1 + (ntot_sml - 1) / NH;
# endif

            // gets memory for gder.
#if !HOST
            cudaMallocManaged(&gder_bb_vec, NB_BIG * nbins_gder * sizeof(float));
            cudaMallocManaged(&gder_ss_vec, NB_SML * nbins_gder * sizeof(float));
            cudaMallocManaged(&gder_bs_vec, NB_BIG * nbins_gder * sizeof(float));
            cudaMallocManaged(&gder_sb_vec, NB_SML * nbins_gder * sizeof(float));
#else
            gder_bb_vec = (float *) malloc(NB_BIG * nbins_gder * sizeof(float));
            gder_ss_vec = (float *) malloc(NB_SML * nbins_gder * sizeof(float));
            gder_bs_vec = (float *) malloc(NB_BIG * nbins_gder * sizeof(float));
            gder_sb_vec = (float *) malloc(NB_SML * nbins_gder * sizeof(float));

#endif
        }

        // positions

        for (int mm = 0; mm < ngrain_tot; mm++) {
            int nn;
            fgets(renglon, sizeof(renglon), fp_in);
            sscanf(renglon, "%d %f %f %f", &nn, &(rr.x), &(rr.y),
                   &(rr.z)); // se debe eliminar esta parte cuando la init config sea correcta

            if (nn != mm) {
                printf("error: mm %d  nn %d no match\n", mm, nn);
                exit(1);
            }
            if (mm < ngrain_big) {
                rr_big_vec[mm] = rr;
                rr_big_raw_vec[mm] = rr;
            } else {
                rr_sml_vec[mm - ngrain_big] = rr;
                rr_sml_raw_vec[mm - ngrain_big] = rr;
            }
        }
        fclose(fp_in);


        // other parameters

        pars.not_finished = 1;
        niter = (int) (run_time / dt);
        ngap = niter / nsamples;
        if (ngap < 1) {
            printf("error: ngap %d\n", ngap);
            exit(1);
        }
        ntrans = (int) (trans_time / dt);
        fprintf(fp_bitac, "niter, ntrans, nsamples, ngap %d %d %d %d\n", niter, ntrans,
                nsamples, ngap);

        //==================================================================================
        // transient and run
        //==================================================================================

        // clear cell vectors and do locate

#if HOST
        set_vec_int_hst(nocup_big_vec, ncell_big3, 0);
        set_vec_int_hst(cell_big_vec, ntot_big, -1);
        cell_locate_hst('b', rr_big_vec, nocup_big_vec, cell_big_vec, pars);
        set_vec_int_hst(nocup_sml_vec, ncell_sml3, 0);
        set_vec_int_hst(cell_sml_vec, ntot_sml, -1);
        cell_locate_hst('s', rr_sml_vec, nocup_sml_vec, cell_sml_vec, pars);
# else
        set_vec_int_dev<<<NB_CELL_BIG3, NH>>>(nocup_big_vec, ncell_big3, 0);
        set_vec_int_dev<<<NB_NTOT_BIG, NH>>>(cell_big_vec, ntot_big, -1);
        cell_locate_dev<<<NB_BIG, NH>>>('b', rr_big_vec, nocup_big_vec, cell_big_vec, pars);
        set_vec_int_dev<<<NB_CELL_SML3, NH>>>(nocup_sml_vec, ncell_sml3, 0);
        set_vec_int_dev<<<NB_NTOT_SML, NH>>>(cell_sml_vec, ntot_sml, -1);
        cell_locate_dev<<<NB_SML, NH>>>('s', rr_sml_vec, nocup_sml_vec, cell_sml_vec, pars);
# endif

        // clear forces

# if HOST
        set_vec_float3_hst(ff_big_vec, ngrain_big, zero);
        set_vec_float_hst(vir_big_vec, ngrain_big, 0.0);
        set_vec_float_hst(pot_big_vec, ngrain_big, 0.0);
        set_vec_float3_hst(ff_sml_vec, ngrain_sml, zero);
        set_vec_float_hst(vir_sml_vec, ngrain_sml, 0.0);
        set_vec_float_hst(pot_sml_vec, ngrain_sml, 0.0);
# else
        set_vec_float3_dev<<<NB_BIG, NH>>>(ff_big_vec, ngrain_big, zero);
        set_vec_float_dev<<<NB_BIG, NH>>>(vir_big_vec, ngrain_big, 0.0);
        set_vec_float_dev<<<NB_BIG, NH>>>(pot_big_vec, ngrain_big, 0.0);
        set_vec_float3_dev<<<NB_SML, NH>>>(ff_sml_vec, ngrain_sml, zero);
        set_vec_float_dev<<<NB_SML, NH>>>(vir_sml_vec, ngrain_sml, 0.0);
        set_vec_float_dev<<<NB_SML, NH>>>(pot_sml_vec, ngrain_sml, 0.0);
# endif

        // forces

# if HOST
        get_forces_same_hst('b', rr_big_vec, ff_big_vec, vir_big_vec, pot_big_vec,
                            nocup_big_vec, cell_big_vec, pars);
        get_forces_same_hst('s', rr_sml_vec, ff_sml_vec, vir_sml_vec, pot_sml_vec,
                            nocup_sml_vec, cell_sml_vec, pars);
        get_forces_diff_hst('b', 's', rr_big_vec, ff_big_vec, rr_sml_vec, vir_big_vec,
                            pot_big_vec, nocup_sml_vec, cell_sml_vec, pars);
        get_forces_diff_hst('s', 'b', rr_sml_vec, ff_sml_vec, rr_big_vec, vir_sml_vec,
                            pot_sml_vec, nocup_big_vec, cell_big_vec, pars);
# else
        get_forces_same_dev<<<NB_BIG, NH>>>('b', rr_big_vec, ff_big_vec, vir_big_vec,
                                            pot_big_vec, nocup_big_vec, cell_big_vec, pars);
        get_forces_same_dev<<<NB_SML, NH>>>('s', rr_sml_vec, ff_sml_vec, vir_sml_vec,
                                            pot_sml_vec, nocup_sml_vec, cell_sml_vec, pars);
        get_forces_diff_dev<<<NB_BIG, NH>>>('b', 's', rr_big_vec, ff_big_vec, rr_sml_vec,
                                            vir_big_vec, pot_big_vec, nocup_sml_vec, cell_sml_vec, pars);
        get_forces_diff_dev<<<NB_SML, NH>>>('s', 'b', rr_sml_vec, ff_sml_vec, rr_big_vec,
                                            vir_sml_vec, pot_sml_vec, nocup_big_vec, cell_big_vec, pars);
#endif

        // initialize Random Number Generators
#if HOST
        gsl_rng *rand = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rand, idum);
#else
        curandState *devStates_big;
        curandState *devStates_sml;
        cudaMalloc(&devStates_big, (size_t) ngrain_big * sizeof(curandState));
        cudaMalloc(&devStates_sml, (size_t) ngrain_sml * sizeof(curandState));

        setup_rng_kernel<<<NB_BIG, NH>>>(devStates_big, idum, ngrain_big);
        setup_rng_kernel<<<NB_SML, NH>>>(devStates_sml, idum + 1, ngrain_sml);
#endif

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////// run ///////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef LOG_SAMPLE
        int sampling_counter = 0;
#endif

        for (ni = -ntrans; ni < niter; ni++) {
            cudaDeviceSynchronize();

            if (ni == 0) {
                cudaMemcpy(rr_big_raw_vec, rr_big_vec, ngrain_big * sizeof(float3),
                           cudaMemcpyHostToHost);
                cudaMemcpy(rr_big_ini_vec, rr_big_vec, ngrain_big * sizeof(float3),
                           cudaMemcpyHostToHost);
                cudaMemcpy(rr_sml_raw_vec, rr_sml_vec, ngrain_sml * sizeof(float3),
                           cudaMemcpyHostToHost);
                cudaMemcpy(rr_sml_ini_vec, rr_sml_vec, ngrain_sml * sizeof(float3),
                           cudaMemcpyHostToHost);
            }

            //Evolution Ermak-McCammon
#if HOST
            update_ermak_hst('b', rr_big_vec, rr_big_raw_vec, ff_big_vec, rand, pars);
            update_ermak_hst('s', rr_sml_vec, rr_sml_raw_vec, ff_sml_vec, rand, pars);
#else
            update_ermak_dev<<<NB_BIG, NH>>>('b', rr_big_vec, rr_big_raw_vec, ff_big_vec, devStates_big, pars);
            update_ermak_dev<<<NB_SML, NH>>>('s', rr_sml_vec, rr_sml_raw_vec, ff_sml_vec, devStates_sml, pars);

            cudaDeviceSynchronize();
#endif

            //calculate time and temperature
            time = dt * (1.0 + (float) ni);
            if (strcmp(temperature_protocol, "linear") == 0)
                pars.temp_set = calculate_temp_linear(T0, Tf, t0, tf, time);
            if (strcmp(temperature_protocol, "sine") == 0)
                pars.temp_set = calculate_temp_sine(T0, Tf, t0, tf, period, time);

            // clear cell vectors and do cell locate

# if HOST
            set_vec_int_hst(nocup_big_vec, ncell_big3, 0);
            set_vec_int_hst(cell_big_vec, ntot_big, -1);
            cell_locate_hst('b', rr_big_vec, nocup_big_vec, cell_big_vec, pars);
            set_vec_int_hst(nocup_sml_vec, ncell_sml3, 0);
            set_vec_int_hst(cell_sml_vec, ntot_sml, -1);
            cell_locate_hst('s', rr_sml_vec, nocup_sml_vec, cell_sml_vec, pars);
# else
            set_vec_int_dev<<<NB_CELL_BIG3, NH>>>(nocup_big_vec, ncell_big3, 0);
            set_vec_int_dev<<<NB_NTOT_BIG, NH>>>(cell_big_vec, ntot_big, -1);
            cell_locate_dev<<<NB_BIG, NH>>>('b', rr_big_vec, nocup_big_vec,
                                            cell_big_vec, pars);
            set_vec_int_dev<<<NB_CELL_SML3, NH>>>(nocup_sml_vec, ncell_sml3, 0);
            set_vec_int_dev<<<NB_NTOT_SML, NH>>>(cell_sml_vec, ntot_sml, -1);
            cell_locate_dev<<<NB_SML, NH>>>('s', rr_sml_vec, nocup_sml_vec,
                                            cell_sml_vec, pars);
# endif

            // clear forces

# if HOST
            set_vec_float3_hst(ff_big_vec, ngrain_big, zero);
            set_vec_float_hst(vir_big_vec, ngrain_big, 0.0);
            set_vec_float_hst(pot_big_vec, ngrain_big, 0.0);
            set_vec_float3_hst(ff_sml_vec, ngrain_sml, zero);
            set_vec_float_hst(vir_sml_vec, ngrain_sml, 0.0);
            set_vec_float_hst(pot_sml_vec, ngrain_sml, 0.0);
# else
            set_vec_float3_dev<<<NB_BIG, NH>>>(ff_big_vec, ngrain_big, zero);
            set_vec_float_dev<<<NB_BIG, NH>>>(vir_big_vec, ngrain_big, 0.0);
            set_vec_float_dev<<<NB_BIG, NH>>>(pot_big_vec, ngrain_big, 0.0);
            set_vec_float3_dev<<<NB_SML, NH>>>(ff_sml_vec, ngrain_sml, zero);
            set_vec_float_dev<<<NB_SML, NH>>>(vir_sml_vec, ngrain_sml, 0.0);
            set_vec_float_dev<<<NB_SML, NH>>>(pot_sml_vec, ngrain_sml, 0.0);
# endif

            // get new forces

# if HOST
            get_forces_same_hst('b', rr_big_vec, ff_big_vec, vir_big_vec, pot_big_vec,
                                nocup_big_vec, cell_big_vec, pars);
            get_forces_same_hst('s', rr_sml_vec, ff_sml_vec, vir_sml_vec, pot_sml_vec,
                                nocup_sml_vec, cell_sml_vec, pars);
            get_forces_diff_hst('b', 's', rr_big_vec, ff_big_vec, rr_sml_vec, vir_big_vec,
                                pot_big_vec, nocup_sml_vec, cell_sml_vec, pars);
            get_forces_diff_hst('s', 'b', rr_sml_vec, ff_sml_vec, rr_big_vec, vir_sml_vec,
                                pot_sml_vec, nocup_big_vec, cell_big_vec, pars);
# else
            get_forces_same_dev<<<NB_BIG, NH>>>('b', rr_big_vec, ff_big_vec,
                                                vir_big_vec, pot_big_vec, nocup_big_vec, cell_big_vec, pars);
            get_forces_same_dev<<<NB_SML, NH>>>('s', rr_sml_vec, ff_sml_vec,
                                                vir_sml_vec, pot_sml_vec, nocup_sml_vec, cell_sml_vec, pars);
            get_forces_diff_dev<<<NB_BIG, NH>>>('b', 's', rr_big_vec, ff_big_vec,
                                                rr_sml_vec, vir_big_vec, pot_big_vec, nocup_sml_vec, cell_sml_vec,
                                                pars);
            get_forces_diff_dev<<<NB_SML, NH>>>('s', 'b', rr_sml_vec, ff_sml_vec,
                                                rr_big_vec, vir_sml_vec, pot_sml_vec, nocup_big_vec, cell_big_vec,
                                                pars);
            cudaDeviceSynchronize();
# endif

            // processing

            if (ni < 0 && ni % ngap == 0) {
                cudaDeviceSynchronize();

                // get nocup_max

                nocup_big_max = nocup_sml_max = 0;
                for (ii = 0; ii < ncell_big3; ii++) {
                    nocup = nocup_big_vec[ii];
                    if (nocup_big_max < nocup) nocup_big_max = nocup;
                }
                for (ii = 0; ii < ncell_sml3; ii++) {
                    nocup = nocup_sml_vec[ii];
                    if (nocup_sml_max < nocup) nocup_sml_max = nocup;
                }
                printf("run ni %d/%d  %.2f%%  --  print %d/%d  --  ocup b s  %d %d\n",
                       ni, niter, 100.0 * ((float) (ni)) / ((float) niter), counter, nsamples,
                       nocup_big_max, nocup_sml_max);
                fflush(stdout);
            }

            // processing
#ifdef LOG_SAMPLE
            if (ni == sampling_idx[sampling_counter]) {
                should_sample = true;
                sampling_counter++;
            }
#else
            if (ni >= 0 && ni % ngap == 0)
                should_sample = true;
#endif

            if (should_sample) {
                cudaDeviceSynchronize();
                counter++;

                // energy sampling
                ene_pot_big = ene_pot_sml = 0.0;
                for (int mm = 0; mm < ngrain_big; mm++) ene_pot_big += 0.5 * pot_big_vec[mm];
                ene_pot_big /= xngrain_big;
                for (int mm = 0; mm < ngrain_sml; mm++) ene_pot_sml += 0.5 * pot_sml_vec[mm];
                ene_pot_sml /= xngrain_sml;

                energy_us[energy_counter] += ene_pot_sml;
                energy_ub[energy_counter] += ene_pot_big;
                time_energy[energy_counter] = time;
                energy_temp[energy_counter] = pars.temp_set;

                energy_counter++;

                // occupation

                nocup_big_max = nocup_sml_max = 0;
                for (ii = 0; ii < ncell_big3; ii++) {
                    nocup = nocup_big_vec[ii];
                    if (nocup_big_max < nocup) nocup_big_max = nocup;
                }
                for (ii = 0; ii < ncell_sml3; ii++) {
                    nocup = nocup_sml_vec[ii];
                    if (nocup_sml_max < nocup) nocup_sml_max = nocup;
                }
                printf("run ni %d/%d  %.2f%%  --  print %d/%d  --  ocup b s  %d %d -- temp %.3f\n",
                       ni, niter, 100.0 * ((float) (ni)) / ((float) niter), counter, nsamples,
                       nocup_big_max, nocup_sml_max, pars.temp_set);
                fflush(stdout);
                if (nocup_big_max >= ntags_big || nocup_sml_max >= ntags_sml) {
                    printf("ntags reached limit\n");
                    fflush(stdout);
                    exit(1);
                }

                // writing snapshots

#if GRO_FLAG
                bool high_render = true;
                char color[20];
                int blue_count, green_count, red_count;
                sprintf(snapshots_fn, "snapshot_%d.pov", counter);
                fp_snaps = fopen(snapshots_fn, "w");
                fprintf(fp_snaps, "#include \"colors.inc\"\n");
                if (high_render) {
                    fprintf(fp_snaps, "global_settings {\n");
                    fprintf(fp_snaps, "\tradiosity{\n");
                    fprintf(fp_snaps, "\t\tpretrace_start 0.08\n");
                    fprintf(fp_snaps, "\t\tpretrace_end   0.01\n");
                    fprintf(fp_snaps, "\t\tcount 150\n");
                    fprintf(fp_snaps, "\t\tnearest_count 10\n");
                    fprintf(fp_snaps, "\t\terror_bound 0.35\n");
                    fprintf(fp_snaps, "\t\trecursion_limit 2\n");
                    fprintf(fp_snaps, "\t\tlow_error_factor 0.5\n");
                    fprintf(fp_snaps, "\t\tgray_threshold 0.0\n");
                    fprintf(fp_snaps, "\t\tminimum_reuse 0.005\n");
                    fprintf(fp_snaps, "\t\tmaximum_reuse 0.2\n");
                    fprintf(fp_snaps, "\t\tbrightness 1\n");
                    fprintf(fp_snaps, "\t\tadc_bailout 0.005\n");
                    fprintf(fp_snaps, "\t}\n}\n");
                }
                fprintf(fp_snaps, "background {White}\n");
                fprintf(fp_snaps, "camera{\n\tangle 90\n\tlocation <9,9,15>\n\tlook_at <9,9,2>\n}\n");
                fprintf(fp_snaps, "light_source{ <25,25,25> color White}\n");

                blue_count = green_count = red_count = 0;

                //print big particles
                for (int mm = 0; mm < ngrain_big; ++mm) {
                    rr = rr_big_vec[mm];
                    diameter = sigma_big;
                    calculate_color(mm, rr_big_vec, rr_sml_vec, color, blue_count, green_count, red_count, pars);
                    if (1.0 < rr.z && rr.z < 3.0)
                        fprintf(fp_snaps, "sphere{ <%8.3f,%8.3f,%8.3f>, %8.3f pigment {%s} finish { specular 0.5 }}\n",
                                rr.x, rr.y, rr.z, diameter * .5, color);
                }

                //print small particles
                for (int mm = 0; mm < ngrain_sml; mm++) {
                    rr = rr_sml_vec[mm];
                    diameter = sigma_sml;
                    if (1.2 < rr.z && rr.z < 2.8)
                        fprintf(fp_snaps,
                                "sphere{ <%8.3f%8.3f%8.3f>, %8.3f pigment {Yellow} finish { specular 0.5 }}\n", rr.x,
                                rr.y, rr.z, diameter * .5);
                }

                fprintf(fp_colors, "%f  %d  %d  %d\n", time, red_count, green_count, blue_count);
                fclose(fp_snaps);
# endif

                //get statistics for gder

                for (int nb = 0; nb < NB_BIG * nbins_gder; nb++) gder_bb_vec[nb] = 0.0f;
                for (int nb = 0; nb < NB_SML * nbins_gder; nb++) gder_ss_vec[nb] = 0.0f;
                for (int nb = 0; nb < NB_BIG * nbins_gder; nb++) gder_bs_vec[nb] = 0.0f;
                for (int nb = 0; nb < NB_SML * nbins_gder; nb++) gder_sb_vec[nb] = 0.0f;

# if HOST
                get_gder_hst('b', 'b', rr_big_vec, rr_big_vec, gder_bb_vec, pars);
                get_gder_hst('b', 's', rr_big_vec, rr_sml_vec, gder_bs_vec, pars);
                get_gder_hst('s', 'b', rr_sml_vec, rr_big_vec, gder_sb_vec, pars);
                get_gder_hst('s', 's', rr_sml_vec, rr_sml_vec, gder_ss_vec, pars);
# else
                get_gder_dev<<<NB_BIG, NH>>>('b', 'b', rr_big_vec, rr_big_vec, gder_bb_vec,
                                             pars);
                get_gder_dev<<<NB_BIG, NH>>>('b', 's', rr_big_vec, rr_sml_vec, gder_bs_vec,
                                             pars);
                get_gder_dev<<<NB_SML, NH>>>('s', 'b', rr_sml_vec, rr_big_vec, gder_sb_vec,
                                             pars);
                get_gder_dev<<<NB_SML, NH>>>('s', 's', rr_sml_vec, rr_sml_vec, gder_ss_vec,
                                             pars);

                cudaDeviceSynchronize();
# endif

                for (jj = 0; jj < nbins_gder; jj++) {
                    for (int nb = 1; nb < NB_BIG; nb++)
                        gder_bb_vec[jj] += gder_bb_vec[jj + nb * nbins_gder];
                    for (int nb = 1; nb < NB_SML; nb++)
                        gder_ss_vec[jj] +=
                                gder_ss_vec[jj + nb * nbins_gder];
                    for (int nb = 1; nb < NB_SML; nb++)
                        gder_sb_vec[jj] +=
                                gder_sb_vec[jj + nb * nbins_gder];
                    for (int nb = 1; nb < NB_BIG; nb++)
                        gder_bs_vec[jj] +=
                                gder_bs_vec[jj + nb * nbins_gder];
                }

                // add to gders arrays to maintain statistics

                if (counter <= 0 || counter > nsamples) {
                    fprintf(stderr, "ERROR: counter=%d fuera de rango\n", counter);
                    exit(EXIT_FAILURE);
                }

                for (int nb = 0; nb < nbins_gder; nb++) {
                    gders[counter - 1][0][nb] += gder_bb_vec[nb];
                    gders[counter - 1][1][nb] += gder_bs_vec[nb];
                    gders[counter - 1][2][nb] += gder_sb_vec[nb];
                    gders[counter - 1][3][nb] += gder_ss_vec[nb];
                }

                // calculate fself

                fself_isotropic(rr_big_raw_vec, rr_big_ini_vec, ngrain_big, fself_big, counter, qmax);
                fself_isotropic(rr_sml_raw_vec, rr_sml_ini_vec, ngrain_sml, fself_sml, counter, qmax);

                // calculate MSD

                msd_big = 0.0;
                for (int mm = 0; mm < ngrain_big; mm++) {
                    drr.x = rr_big_raw_vec[mm].x - rr_big_ini_vec[mm].x;
                    drr.y = rr_big_raw_vec[mm].y - rr_big_ini_vec[mm].y;
                    drr.z = rr_big_raw_vec[mm].z - rr_big_ini_vec[mm].z;
                    msd_big += drr.x * drr.x + drr.y * drr.y + drr.z * drr.z;
                }
                msd_big /= xngrain_big;

                msd_sml = 0.0;
                for (int mm = 0; mm < ngrain_sml; mm++) {
                    drr.x = rr_sml_raw_vec[mm].x - rr_sml_ini_vec[mm].x;
                    drr.y = rr_sml_raw_vec[mm].y - rr_sml_ini_vec[mm].y;
                    drr.z = rr_sml_raw_vec[mm].z - rr_sml_ini_vec[mm].z;
                    msd_sml += drr.x * drr.x + drr.y * drr.y + drr.z * drr.z;
                }
                msd_sml /= xngrain_sml;


                time_msd[counter - 1] = time;
                msd_b[counter - 1] = msd_big;
                msd_s[counter - 1] = msd_sml;

                // calcula presion

                virial = 0.0;
                for (int mm = 0; mm < ngrain_big; mm++) virial += vir_big_vec[mm];
                for (int mm = 0; mm < ngrain_sml; mm++) virial += vir_sml_vec[mm];
                big_z = 1.0 + virial / (3.0 * xngrain_tot * pars.temp_set); // compressibility 3D
                pressure[counter - 1] += big_z;
                time_pressure[counter - 1] = time;

                // print configs for other decade
#ifdef LOG_SAMPLE
                if ((counter - 1) % points_per_decade == 1)
                    if (print_decades_configs && (number_of_tw <= number_of_tws)) {
                        sprintf(infile, "../configs/decade%d/init_config_%d", number_of_tw, i_config);
                        fp_out = fopen(infile, "w");

                        fprintf(fp_out, "%f %f    phi b, s\n", phi_big, phi_sml);
                        fprintf(fp_out, "%f    side\n", side);
                        fprintf(fp_out, "%d %d %d    ngr b, s, tot\n", ngrain_big, ngrain_sml, ngrain_tot);
                        fprintf(fp_out, "%f %f    sigma b, s\n", sigma_big, sigma_sml);
                        fprintf(fp_out, "%f %f    mass b, s\n", mass_big, mass_sml);
                        fprintf(fp_out, "\n");

                        for (int mm = 0; mm < ngrain_big; mm++) {
                            rr = rr_big_vec[mm];
                            fprintf(fp_out, "%d %f %f %f\n", mm, rr.x, rr.y, rr.z);
                        }
                        for (int mm = 0; mm < ngrain_sml; mm++) {
                            rr = rr_sml_vec[mm];
                            fprintf(fp_out, "%d %f %f %f\n", (ngrain_big + mm),
                                    rr.x, rr.y, rr.z);
                        }
                        fclose(fp_out);

                        number_of_tw++;
                    }
#endif

                // restart condition
                should_sample = false;
            }
        }
        cudaDeviceSynchronize();

        printf("Finished step %d of %d\n", i_config, n_configs);
    }


    //normaliza gder

    for (int i = 0; i < nsamples; ++i)
        for (int nb = 0; nb < nbins_gder; nb++) {
            xnb = (float) nb;

            shell_vol = 4.0 * PI * bin_size_gder * bin_size_gder * bin_size_gder *
                        ((1.0 / 3.0) + xnb * xnb + xnb);

            sigma = sigma_big;
            vol_free = volume - (4.0 / 3.0) * PI * sigma * sigma * sigma;
            gders[i][0][nb] = gders[i][0][nb] / (shell_vol);
            gders[i][0][nb] /= (xngrain_big * (xngrain_big - 1.0) / vol_free);

            sigma = 0.5 * (sigma_big + sigma_sml);
            vol_free = volume - (4.0 / 3.0) * PI * sigma * sigma * sigma;
            gders[i][1][nb] = gders[i][1][nb] / (shell_vol);
            gders[i][2][nb] = gders[i][2][nb] / (shell_vol);
            gders[i][1][nb] /= (xngrain_big * xngrain_sml / vol_free);
            gders[i][2][nb] /= (xngrain_sml * xngrain_big / vol_free);

            sigma = sigma_sml;
            vol_free = volume - (4.0 / 3.0) * PI * sigma * sigma * sigma;
            gders[i][3][nb] = gders[i][3][nb] / (shell_vol);
            gders[i][3][nb] /= (xngrain_sml * (xngrain_sml - 1.0) / vol_free);
        }

    for (int i = 0; i < nsamples; ++i)
        for (int j = 0; j < 4; ++j)
            for (int nb = 0; nb < nbins_gder; ++nb)
                gders[i][j][nb] /= n_configs;

    for (int i = 0; i < nsamples; ++i) {
        energy_ub[i] /= n_configs;
        energy_us[i] /= n_configs;
    }


    //imprime resultados en files
    printf("printing files\n");

    for (int i = 0; i < nsamples; ++i) {
        if (strcmp(sample_type, "log") == 0)
            sprintf(gder_fn, "results/decade%d/gder_%d.out", current_decade, i);
        else
            sprintf(gder_fn, "results/gder_%d.out", i);

        fp_gder = fopen(gder_fn, "w");
        if (fp_gder == nullptr) {
            printf("Verify file path: %s", gder_fn);
            fflush(stdout);
            exit(-2);
        }

        for (int nb = 0; nb < nbins_gder; nb++) {
            dist = (0.5 + (float) nb) * bin_size_gder;
            fprintf(fp_gder, "%f  %f  %f  %f  %f\n", dist, gders[i][0][nb], gders[i][1][nb],
                    gders[i][2][nb], gders[i][3][nb]);
        }
        fclose(fp_gder);
    }

    if (strcmp(sample_type, "log") == 0)
        sprintf(energy_fn, "results/decade%d/energies_averaged.out", current_decade);
    else
        sprintf(energy_fn, "results/energies_averaged.out");

    fp_energ = fopen(energy_fn, "w");
    if (fp_energ == nullptr) {
        printf("Verify file path: %s", energy_fn);
        fflush(stdout);
        exit(-2);
    }

    for (int i = 0; i < nsamples; ++i) {
        fprintf(fp_energ, "%f  %f  %f  %f\n", time_energy[i], energy_ub[i],
                energy_us[i], energy_temp[i]);
    }
    fclose(fp_energ);

    for (int i = 0; i < nsamples; ++i) {
        fself_big[i] /= n_configs;
        fself_sml[i] /= n_configs;
    }

    if (strcmp(sample_type, "log") == 0)
        sprintf(fself_fn, "results/decade%d/fself.out", current_decade);
    else
        sprintf(fself_fn, "results/fself.out");

    fp_fself = fopen(fself_fn, "w");
    if (fp_fself == nullptr) {
        printf("Verify file path: %s", fself_fn);
        fflush(stdout);
        exit(-2);
    }
    for (int i = 0; i < nsamples; ++i)
        fprintf(fp_fself, "%f  %f  %f\n", time_msd[i], fself_big[i], fself_sml[i]);

    fclose(fp_fself);

    if (strcmp(sample_type, "log") == 0)
        sprintf(msd_fn, "results/decade%d/msd_averaged.out", current_decade);
    else
        sprintf(msd_fn, "results/msd_averaged.out");

    fp_msd = fopen(msd_fn, "w");
    if (fp_msd == nullptr) {
        printf("Verify file path: %s", msd_fn);
        fflush(stdout);
        exit(-2);
    }

    for (int i = 0; i < nsamples; ++i) {
        fprintf(fp_msd, "%f  %f  %f\n", time_msd[i], msd_b[i], msd_s[i]);
    }
    fclose(fp_msd);

    // pressure
    if (strcmp(sample_type, "log") == 0)
        sprintf(press_fn, "results/decade%d/pressure.out", current_decade);
    else
        sprintf(press_fn, "results/pressure.out");

    fp_press = fopen(press_fn, "w");
    if (fp_press == nullptr) {
        printf("Verify file path: %s", press_fn);
        fflush(stdout);
        exit(-2);
    }

    for (int i = 0; i < nsamples; ++i) {
        pressure[i] /= n_configs;
        fprintf(fp_press, "%f  %f\n", time_pressure[i], pressure[i]);
    }
    fclose(fp_press);


    // close files, release stuff and finish

#if GRO_FLAG
    fclose(fp_colors);
#endif

    fclose(fp_bitac);

    cudaFree(rr_big_vec);
    cudaFree(rr_big_raw_vec);
    cudaFree(rr_big_ini_vec);
    cudaFree(ff_big_vec);
    cudaFree(vir_big_vec);
    cudaFree(pot_big_vec);

    cudaFree(nocup_big_vec);
    cudaFree(cell_big_vec);

    cudaFree(rr_sml_vec);
    cudaFree(rr_sml_raw_vec);
    cudaFree(rr_sml_ini_vec);
    cudaFree(ff_sml_vec);
    cudaFree(vir_sml_vec);
    cudaFree(pot_sml_vec);

    cudaFree(nocup_sml_vec);
    cudaFree(cell_sml_vec);

    cudaFree(gder_bb_vec);
    cudaFree(gder_bs_vec);
    cudaFree(gder_sb_vec);
    cudaFree(gder_ss_vec);


    for (int i = 0; i < nsamples; ++i) {
        for (int j = 0; j < 4; ++j) free(gders[i][j]);
        free(gders[i]);
    }
    free(gders);

    free(fself_big);
    free(fself_sml);

    free(energy_us);
    free(energy_ub);
    free(msd_b);
    free(msd_s);

    free(energy_temp);
    free(time_energy);
    free(time_msd);


    printf("TERMINADO\n");
    exit(0);
}
