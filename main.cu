/* Elaborado por Marco Antonio Ramirez Guizar
Codigo en C++ paralelizado con nvidia-CUDA para simular un sistema bidisperso de particulas
para estudiar las fuerzas efectivas bajo el esquema de contraccion de fuerzas en un sistema
en enfriamiento o calentamiento.
*/
# include "encabezados.h"
#include "polydispersity_IO.h"

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

    PolydispersityParameters poly{};
    SystemParameters pars{};
    float3 rr, drr, zero;
    double run_time, trans_time, dt;
    float side, cutoff, bin_size_gder, qmax, xngrain_big,
            shell_vol, vol_free, msd_big,
            msd_sml, dist, aux, sigma, sigma_big, sigma_sml, mass_big, mass_sml,
            range_gder, xngrain_tot, big_z,
            virial, xnb, T0, Tf, t0, tf, time, period,
            volume;
    float ***gders, *energy_temp,
            *time_energy, *time_msd;
    int ngrain_tot, number_species, ni, niter, ntrans, ngap, idum, ntags_big, ntags_sml,
            nsamples, ii, jj, nbins_gder, counter, ngap_rescaling, n_configs;
    int NH;
    double minimal_decade;
    bool should_sample, print_decades_configs;
    int number_of_tws, number_of_tw, current_decade;

# if !HOST
    int NB_CELL_BIG3, NB_CELL_SML3, NB_NTOT_BIG, NB_NTOT_SML;
# endif
    char renglon[200], infile[80], gder_fn[80], energy_fn[80], msd_fn[80], snapshots_fn[80],
            press_fn[80], fself_fn[80];
    char temperature_protocol[40], sample_type[40], distribution[40];
    FILE *fp_bitac, *fp_snaps, *fp_energ, *fp_gder, *fp_msd, *fp_in, *fp_out,
            *fp_press, *fp_data, *fp_colors, *fp_fself;
    fp_data = fopen("wca_aging_log.data", "r");
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
    pars.temperature = T0;
    pars.nbins_gder = nbins_gder;
    pars.range_gder = range_gder;

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
    number_species = 1;
    //Values to calculate size of arrays
    niter = (int) (run_time / dt);
    ngap = niter / nsamples;

    //g_ij(r) bidisperse
    gders = (float ***) malloc(nsamples * sizeof(float **));
    for (int i = 0; i < nsamples; ++i) {
        gders[i] = (float **) malloc(1 * sizeof(float *));
        for (int j = 0; j < 1; ++j) {
            gders[i][j] = (float *) malloc(nbins_gder * sizeof(float));
            for (int k = 0; k < nbins_gder; ++k) gders[i][j][k] = 0.0;
        }
    }

    //fself
    auto *fself_big = (float *) calloc(nsamples, sizeof(float));

    //energy and msd and pressure
    time_energy = (float *) malloc(nsamples * sizeof(float));

    float **potential_energy;
    potential_energy = (float ** ) malloc(number_species * sizeof(float *));
    for (int i = 0; i < number_species; ++i)
        potential_energy[i] = (float *) calloc(nsamples, sizeof(float));

    energy_temp = (float *) calloc(nsamples, sizeof(float));

    time_msd = (float *) malloc(nsamples * sizeof(float));

    float **msd;
    msd = (float **) malloc(number_species * sizeof(float *));
    for (int i = 0; i < number_species; ++i)
        msd[i] = (float *) calloc(nsamples, sizeof(float));

    auto time_pressure = (float *) malloc(nsamples * sizeof(float));
    auto pressure = (float *) calloc(nsamples, sizeof(float));

#if GRO_FLAG
    fp_colors = fopen("color_statistics.out", "w");
#endif
    ////////////////////////////Starts loop//////////////////////////////////////////////////////////////////////
    xngrain_big = 0.0f;
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
            exit(-2);
        }
        scan_parameters_init_config(fp_in, poly, pars, distribution);

        // almacena datos iniciales
        number_species = pars.number_species;
        ngrain_tot = pars.ngrain_total;
        side = pars.side;

            // Memory allocation

    Grains *grain;
    cudaMallocManaged(&grain, number_species * sizeof(Grains));

    auto NB = (int *) malloc(number_species * sizeof(int));

    for (int i = 0; i < number_species; ++i) {

        grain[i].number_particles = pars.N[i];

        cudaMallocManaged(&(grain[i].rr), grain[i].number_particles * sizeof(float3));
        cudaMallocManaged(&(grain[i].vv), grain[i].number_particles * sizeof(float3));
        cudaMallocManaged(&(grain[i].ff), grain[i].number_particles * sizeof(float3));
        cudaMallocManaged(&(grain[i].diameter), grain[i].number_particles * sizeof(float));
        cudaMallocManaged(&(grain[i].mass), grain[i].number_particles * sizeof(float));

        cudaMallocManaged(&(grain[i].virial), grain[i].number_particles * sizeof(float));
        cudaMallocManaged(&(grain[i].potential), grain[i].number_particles * sizeof(float));

        // for msd
        cudaMallocManaged(&(grain[i].rr0), grain[i].number_particles * sizeof(float3));
        cudaMallocManaged(&(grain[i].rr_raw), grain[i].number_particles * sizeof(float3));

    }

    // bloques e hilos

    NH = 256;
    pars.NH = NH;

    //number of blocks necessary to include all ith particles
    for (int i = 0; i < number_species; ++i) NB[i] = 1 + (grain[i].number_particles - 1) / NH;

    // Radial distribution functions g_ij(r)


    float ***gr;

    cudaMallocManaged(&gr, number_species * sizeof(float **));

    for (int i = 0; i < number_species; ++i) {
        cudaMallocManaged(&(gr[i]), number_species * sizeof(float *));
        for (int j = 0; j < number_species; ++j)
            cudaMallocManaged(&(gr[i][j]), (NB[i] * nbins_gder) * sizeof(float));
    }


    // read vectors for grains
    int nparticles = 0;
    for (int i = 0; i < number_species; ++i) {

        for (int mm = nparticles; mm < nparticles + pars.N[i]; ++mm) {

            int nn;
            float3 rr, vv;
            float diameter, mass;
            fgets(renglon, sizeof(renglon), fp_in);
            sscanf(renglon, "%d %f %f %f %f %f %f %f %f", &nn, &(rr.x), &(rr.y), &(rr.z),
                   &(vv.x), &(vv.y), &(vv.z), &diameter, &mass);
            if (nn != mm) {
                printf("error: %d != %d\n", nn, mm);
                exit(-1);
            }
            grain[i].rr[mm - nparticles] = rr;
            grain[i].vv[mm - nparticles] = vv;
            grain[i].diameter[mm - nparticles] = diameter;
            grain[i].mass[mm - nparticles] = mass;
        }

        nparticles += pars.N[i];
    }

    fclose(fp_in);

    // almacena datos iniciales

    pars.side = side;
    volume = side * side * side;
    xngrain_big += (float) pars.N[0];
    fp_bitac = fopen("bitacora", "w");

    // write in stdin
    print_parameters_simulation(stdin, pars, infile, trans_time, run_time, nsamples);
    print_parameters_init_config(stdin, poly, pars, distribution);

    // write in bitacora
    print_parameters_simulation(fp_bitac, pars, infile, trans_time, run_time, nsamples);
    print_parameters_init_config(fp_bitac, poly, pars, distribution);

    // for gder set it in zero.

    for (int i = 0; i < number_species; ++i)
        for (int j = 0; j < number_species; ++j)
            for (int nbin = 0; nbin < NB[i] * nbins_gder; ++nbin)
                gr[i][j][nbin] = 0.0f;


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


        // clear forces

# if HOST

        for (int i = 0; i < number_species; ++i) {
            set_vec_float3_hst(grain[i].ff, grain[i].number_particles, zero);
            set_vec_float_hst(grain[i].virial, grain[i].number_particles, 0.0);
            set_vec_float_hst(grain[i].potential, grain[i].number_particles, 0.0);
        }

# else
        for (int i = 0; i < number_species; ++i) {
            set_vec_float3_dev<<<NB[i], NH>>>(grain[i].ff, grain[i].number_particles, zero);
            set_vec_float_dev<<<NB[i], NH>>>(grain[i].virial, grain[i].number_particles, 0.0);
            set_vec_float_dev<<<NB[i], NH>>>(grain[i].potential, grain[i].number_particles, 0.0);
        }
# endif

        // forces

# if HOST

        for (int i = 0; i < number_species; ++i) {

            get_forces_same_hst(grain[i], pars);

            for (int j = 0; j < number_species; ++j) {
                if (i == j) continue;
                get_forces_diff_hst(grain[i], grain[j], pars);
            }
        }

# else

        for (int i = 0; i < number_species; ++i) {

            get_forces_same_dev<<<NB[i], NH>>>(grain[i], pars);

            for (int j = 0; j < number_species; ++j) {
                if (i == j) continue;
                get_forces_diff_dev<<<NB[i], NH>>>(grain[i], grain[j], pars);
            }
        }

# endif

        // initialize Random Number Generators
#if HOST
        gsl_rng *rand = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rand, idum);
#else

        curandState **devStates;
        devStates = (curandState **) malloc(sizeof(curandState *));
        for (int i = 0; i < number_species; ++i) {
            cudaMalloc(&devStates[i], (size_t) pars.N[i] * sizeof(curandState));
            setup_rng_kernel<<<NB[i], NH>>>(devStates[i], idum + i, pars.N[i]);
        }

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
                for (int i = 0; i < number_species; ++i) {
                    cudaMemcpy(grain[i].rr0, grain[i].rr, grain[i].number_particles * sizeof(float3), cudaMemcpyHostToHost);
                    cudaMemcpy(grain[i].rr_raw, grain[i].rr, grain[i].number_particles * sizeof(float3),
                               cudaMemcpyHostToHost);
                }
            }

            //Evolution Ermak-McCammon
#if HOST
            for (int i = 0; i < number_species; ++i)
                update_ermak_hst(grain[i], pars, rand);

#else
            for (int i = 0; i < number_species; ++i)
                update_ermak_dev<<<NB[i], NH>>>(grain[i], pars, devStates[i]);

            cudaDeviceSynchronize();
#endif

            //calculate time and temperature
            time = dt * (1.0 + (float) ni);
            if (strcmp(temperature_protocol, "linear") == 0)
                pars.temperature = calculate_temp_linear(T0, Tf, t0, tf, time);
            if (strcmp(temperature_protocol, "sine") == 0)
                pars.temperature = calculate_temp_sine(T0, Tf, t0, tf, period, time);


            // clear forces

# if HOST

            for (int i = 0; i < number_species; ++i) {
                set_vec_float3_hst(grain[i].ff, grain[i].number_particles, zero);
                set_vec_float_hst(grain[i].virial, grain[i].number_particles, 0.0);
                set_vec_float_hst(grain[i].potential, grain[i].number_particles, 0.0);
            }

# else

            for (int i = 0; i < number_species; ++i) {
                set_vec_float3_dev<<<NB[i], NH>>>(grain[i].ff, grain[i].number_particles, zero);
                set_vec_float_dev<<<NB[i], NH>>>(grain[i].virial, grain[i].number_particles, 0.0);
                set_vec_float_dev<<<NB[i], NH>>>(grain[i].potential, grain[i].number_particles, 0.0);
            }

# endif

            // get new forces

# if HOST

            for (int i = 0; i < number_species; ++i) {

                get_forces_same_hst(grain[i], pars);

                for (int j = 0; j < number_species; ++j) {
                    if (i == j) continue;
                    get_forces_diff_hst(grain[i], grain[j], pars);
                }
            }

# else
            for (int i = 0; i < number_species; ++i) {

                get_forces_same_dev<<<NB[i], NH>>>(grain[i], pars);

                for (int j = 0; j < number_species; ++j) {
                    if (i == j) continue;
                    get_forces_diff_dev<<<NB[i], NH>>>(grain[i], grain[j], pars);
                }
            }
            cudaDeviceSynchronize();

# endif

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

                printf("run ni %d/%d  %.2f%%  --  print %d/%d  --\n",
                   ni, niter, 100.0 * ((float) (ni)) / ((float) niter), counter, nsamples);

                // energy sampling
                for (int i = 0; i < number_species; ++i) {
                    float potential = 0;
                    for (int mm = 0; mm < grain[i].number_particles; ++mm)
                        potential += 0.5f * grain[i].potential[mm];

                    potential /= (float) grain[i].number_particles;
                    potential_energy[i][energy_counter] = potential;
                }

                // temperature sampling
                time_energy[energy_counter] = time;
                energy_temp[energy_counter] = pars.temperature;

                energy_counter++;

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

                // reset gder before add new statistics

                for (int i = 0; i < number_species; ++i)
                    for (int j = 0; j < number_species; ++j)
                        for (int nbin = 0; nbin < NB[i] * nbins_gder; ++nbin)
                            gr[i][j][nbin] = 0.0f;


                //get statistics for gder
# if HOST

                for (int i = 0; i < number_species; ++i)
                    for (int j = 0; j < number_species; ++j)
                        get_gder_hst(i, j, grain[i].rr, grain[j].rr, gr[i][j], pars);

# else
                // No funciona en device
                /*for (int i = 0; i < number_species; ++i)
                    for (int j = 0; j < number_species; ++j)
                        get_gder_dev<<<NB[i], NH>>>(i, j, grain[i].rr, grain[j].rr, gr[i][j], pars);*/

                for (int i = 0; i < number_species; ++i)
                    for (int j = 0; j < number_species; ++j)
                        get_gder_hst(i, j, grain[i].rr, grain[j].rr, gr[i][j], pars);


# endif

                for (int i = 0; i < number_species; ++i)
                    for (int j = 0; j < number_species; ++j)
                        for (int nbin = 0; nbin < nbins_gder; ++nbin)
                            for (int number_block = 0; number_block < NB[i]; ++number_block)
                                gr[i][j][nbin] += gr[i][j][nbin + number_block * nbins_gder];

                // add to gders arrays to maintain statistics

                if (counter <= 0 || counter > nsamples) {
                    fprintf(stderr, "ERROR: counter=%d fuera de rango\n", counter);
                    exit(EXIT_FAILURE);
                }

                for (int nb = 0; nb < nbins_gder; nb++)
                    gders[counter - 1][0][nb] += gr[0][0][nb];


                // calculate fself

                fself_isotropic(grain[0].rr_raw, grain[0].rr0, grain[0].number_particles, fself_big, counter, qmax);

                // calculate MSD

                for ( int i = 0; i < number_species; ++i) {
                    float msd_value = 0.0f;
                    for (int mm= 0; mm < grain[i].number_particles; ++mm) {
                        drr.x = grain[i].rr_raw[mm].x - grain[i].rr0[mm].x;
                        drr.y = grain[i].rr_raw[mm].y - grain[i].rr0[mm].y;
                        drr.z = grain[i].rr_raw[mm].z - grain[i].rr0[mm].z;
                        msd_value += drr.x * drr.x + drr.y * drr.y + drr.z * drr.z;
                    }
                    msd_value /= (float) grain[i].number_particles;

                    msd[i][counter - 1] = msd_value;


                }

                time_msd[counter - 1] = time;

                // calcula presion

                virial = 0.0;
                for (int i = 0; i < number_species; ++i)
                    for (int mm = 0; mm < grain[i].number_particles; ++mm)
                        virial += grain[i].virial[mm];

                big_z = 1.0f + virial / (3.0f * xngrain_tot * pars.temperature); // compressibility 3D
                pressure[counter - 1] += big_z;
                time_pressure[counter - 1] = time;

                // print configs for other decade
#ifdef LOG_SAMPLE
                if ((counter - 1) % points_per_decade == 1)
                    if (print_decades_configs && (number_of_tw <= number_of_tws)) {
                        sprintf(infile, "../configs/decade%d/init_config_%d", number_of_tw, i_config);
                        fp_out = fopen(infile, "w");

                        print_parameters_init_config(fp_out, poly, pars, distribution);
                        print_particle_properties(fp_out, grain, pars.number_species);

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
    xngrain_big /= (float) n_configs;

    for (int counter = 0; counter < nsamples; ++counter)
        for (int nb = 0; nb < nbins_gder; nb++) {
            xnb = (float) nb;

            shell_vol = 4.0 * PI * bin_size_gder * bin_size_gder * bin_size_gder *
                        ((1.0 / 3.0) + xnb * xnb + xnb);

            sigma = sigma_big;
            vol_free = volume - (4.0 / 3.0) * PI * sigma * sigma * sigma;
            gders[counter][0][nb] = gders[counter][0][nb] / (shell_vol);
            gders[counter][0][nb] /= (xngrain_big * (xngrain_big - 1.0) / vol_free);

        }

    for (int i = 0; i < nsamples; ++i)
        for (int j = 0; j < 1; ++j)
            for (int nb = 0; nb < nbins_gder; ++nb)
                gders[i][j][nb] /= n_configs;

    for (int counter = 0; counter < nsamples; ++counter) {
        for (int j = 0; j < number_species; ++j)
            potential_energy[j][counter] /= n_configs;

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
            fprintf(fp_gder, "%f  %f\n", dist, gders[i][0][nb]);
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
        fprintf(fp_energ, "%f  %f  %f\n", time_energy[i], potential_energy[0][i], energy_temp[i]);
    }
    fclose(fp_energ);

    for (int i = 0; i < nsamples; ++i)
        fself_big[i] /= n_configs;

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
        fprintf(fp_fself, "%f  %f\n", time_msd[i], fself_big[i]);

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
        fprintf(fp_msd, "%f  %f\n", time_msd[i], msd[0][i]);
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


    for (int i = 0; i < nsamples; ++i) {
        for (int j = 0; j < 1; ++j) free(gders[i][j]);
        free(gders[i]);
    }
    free(gders);

    free(fself_big);
    free(potential_energy);
    free(msd);
    free(pressure);
    free(time_pressure);


    free(energy_temp);
    free(time_energy);
    free(time_msd);


    printf("TERMINADO\n");
    exit(0);
}
