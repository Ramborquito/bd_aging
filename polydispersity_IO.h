//
// Created by marcoramirez on 3/6/25.
//
#pragma once
#ifndef SETCONFIG_BIN_POLYDISPERSITY_IO_H
#define SETCONFIG_BIN_POLYDISPERSITY_IO_H

#endif //SETCONFIG_BIN_POLYDISPERSITY_IO_H


#include "encabezados.h"

void scan_parameters_set_config(FILE *data_file, PolydispersityParameters &poly, char *distribution, SystemParameters &pars,
                                char *init_conf_file) {

    char renglon[200], *ptr;

    printf("Distribution ?\n");
    fgets(renglon, sizeof(renglon), data_file);
    sscanf(renglon, "%s", distribution);

    printf("number of species ?\n");
    fgets(renglon, sizeof(renglon), data_file);
    sscanf(renglon, "%d", &pars.number_species);

    // initialize variables depending of the number of species
    poly.phi = (float *) malloc(pars.number_species * sizeof(float));
    poly.diameter = (float *) malloc(pars.number_species * sizeof(float));
    poly.polydispersity = (float *) malloc(pars.number_species * sizeof(float));
    poly.max_polydispersity = (float *) malloc(pars.number_species * sizeof(float));

    printf("phi_i ?\n");
    fgets(renglon, sizeof(renglon), data_file);
    ptr = renglon;
    for (int i = 0; i < pars.number_species; ++i) poly.phi[i] = strtod(ptr, &ptr);

    printf("diameter_i ?\n");
    fgets(renglon, sizeof(renglon), data_file);
    ptr = renglon;
    for (int i = 0; i < pars.number_species; ++i) poly.diameter[i] = strtod(ptr, &ptr);

    printf("polydispersity: big, sml ?\n");
    fgets(renglon, sizeof(renglon), data_file);
    ptr = renglon;
    for (int i = 0; i < pars.number_species; ++i) poly.polydispersity[i] = strtod(ptr, &ptr);

    printf("max polydispersity: ?\n");
    fgets(renglon, sizeof(renglon), data_file);
    ptr = renglon;
    for (int i = 0; i < pars.number_species; ++i) poly.max_polydispersity[i] = strtod(ptr, &ptr);

    printf("ngrain ?\n");
    fgets(renglon, sizeof(renglon), data_file);
    sscanf(renglon, "%d", &pars.ngrain_total);

    printf("temp ?\n");
    fgets(renglon, sizeof(renglon), data_file);
    sscanf(renglon, "%f", &pars.temperature);

    printf("idum ?\n");
    fgets(renglon, sizeof(renglon), data_file);
    sscanf(renglon, "%d", &pars.idum);

    printf("init_config file ?\n");
    fgets(renglon, sizeof(renglon), data_file);
    sscanf(renglon, "%s", init_conf_file);

}

void print_parameters_init_config(FILE *fp, PolydispersityParameters poly, SystemParameters pars,
                                  char *distribution) {

    fprintf(fp, "%d  ngrain t, b, s\n",
            pars.ngrain_total);

    fprintf(fp, "%f  %f  side, temperature\n", pars.side, pars.temperature);

    fprintf(fp, "%d  %s  number of species, distribution\n", pars.number_species, distribution);

    for (int i = 0; i < pars.number_species; ++i) fprintf(fp, "%f  ", poly.phi[i]);
    fprintf(fp, " phis\n");

    for (int i = 0; i < pars.number_species; ++i) fprintf(fp, "%d  ", pars.N[i]);
    fprintf(fp, " N_i\n");

    for (int i = 0; i < pars.number_species; ++i) fprintf(fp, "%f  ", poly.diameter[i]);
    fprintf(fp, " mean diameters\n");

    for (int i = 0; i < pars.number_species; ++i) fprintf(fp, "%f  ", poly.polydispersity[i]);
    fprintf(fp, " polydispersities\n");

    for (int i = 0; i < pars.number_species; ++i) fprintf(fp, "%f  ", poly.max_polydispersity[i]);
    fprintf(fp, " max_polydispersities\n");

}

void print_particle_properties(FILE *output_file, const Grains *grain_vec, const int number_species) {

    int nparticles = 0;
    for (int i = 0; i < number_species; ++i) {

        for (int mm = nparticles; mm < grain_vec[i].number_particles; mm++) {
            fprintf(output_file, "%d %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n",
                    mm,
                    grain_vec[i].rr[mm].x, grain_vec[i].rr[mm].y, grain_vec[i].rr[mm].z,
                    grain_vec[i].vv[mm].x, grain_vec[i].vv[mm].y, grain_vec[i].vv[mm].z,
                    grain_vec[i].diameter[mm], grain_vec[i].mass[mm]
            );
        }
        nparticles += grain_vec[i].number_particles;
    }
}

void scan_parameters_init_config(FILE *fp, PolydispersityParameters &poly, SystemParameters &pars,
                                 char *distribution) {

    char renglon[200], *ptr;

    fgets(renglon, sizeof(renglon), fp);
    sscanf(renglon, "%d", &pars.ngrain_total);

    fgets(renglon, sizeof(renglon), fp);
    sscanf(renglon, "%f  %f", &pars.side, &pars.temperature);

    fgets(renglon, sizeof(renglon), fp);
    sscanf(renglon, "%d  %s", &pars.number_species, distribution);

    // Initialize memory for polydispersity parameters
    poly.phi = (float*)malloc(pars.number_species * sizeof(float));
    poly.diameter = (float*)malloc(pars.number_species * sizeof(float));
    poly.polydispersity = (float*)malloc(pars.number_species * sizeof(float));
    poly.max_polydispersity = (float*)malloc(pars.number_species * sizeof(float));

    // initialize grains parameters

    pars.N = (int*)malloc(pars.number_species * sizeof(int));

    fgets(renglon, sizeof(renglon), fp);
    ptr = renglon;
    for (int i = 0; i < pars.number_species; ++i) poly.phi[i] = strtod(ptr, &ptr);

    fgets(renglon, sizeof(renglon), fp);
    ptr = renglon;
    for (int i = 0; i < pars.number_species; ++i) pars.N[i] = (int)strtol(ptr, &ptr, 10);

    fgets(renglon, sizeof(renglon), fp);
    ptr = renglon;
    for (int i = 0; i < pars.number_species; ++i) poly.diameter[i] = strtof(ptr, &ptr);

    fgets(renglon, sizeof(renglon), fp);
    ptr = renglon;
    for (int i = 0; i < pars.number_species; ++i) poly.polydispersity[i] = strtof(ptr, &ptr);

    fgets(renglon, sizeof(renglon), fp);
    ptr = renglon;
    for (int i = 0; i < pars.number_species; ++i) poly.max_polydispersity[i] = strtof(ptr, &ptr);

}

//********************* Falta hacer bien
void scan_particle_properties_init_config(FILE *fp, Grains grain_vec, const int ngrain_total) {

    char renglon[200];
    int nn;
    float rx, ry, rz, vx, vy, vz;
    float diameter, mass;

    for (int mm = 0; mm < ngrain_total; mm++)
    {
        fgets(renglon, sizeof(renglon), fp);
        sscanf(renglon, "%d %f %f %f %f %f %f %f %f", &nn, &rx, &ry, &rz,
               &vx, &vy, &vz, &diameter, &mass);
        if (nn != mm)
        {
            printf("error_02: nn %d differs mm %d\n", nn, mm);
            exit (1);
        }
        grain_vec.rr[mm].x = rx;
        grain_vec.rr[mm].y = ry;
        grain_vec.rr[mm].z = rz;
        grain_vec.vv[mm].x = vx;
        grain_vec.vv[mm].y = vy;
        grain_vec.vv[mm].z = vz;
        grain_vec.diameter[mm] = diameter;
        grain_vec.mass[mm] = mass;
    }

}


void scan_parameters_eliminate_overlap(FILE *data_file, char *init_conf_file_in, char *bitac_file,
                                       char *energy_file, char *snapshot_file, char *init_conf_file_out,
                                       OverlapParameters &overlap_parameters,  float &time_run, int &idum){

    char renglon[200];

    printf("init_conf_file_in?\n");
    fgets(renglon, sizeof(renglon), data_file);
    sscanf(renglon, "%s", init_conf_file_in);

    printf("bitac_file, energ_file?\n");
    fgets(renglon, sizeof(renglon), data_file);
    sscanf(renglon, "%s %s", bitac_file, energy_file);

    printf("snapshot_file, init_conf_file_out?\n");
    fgets(renglon, sizeof(renglon), data_file);
    sscanf(renglon, "%s %s", snapshot_file, init_conf_file_out);

    printf("sdp kappa, gamma?\n");
    fgets(renglon, sizeof(renglon), data_file);
    sscanf(renglon, "%f %f", &overlap_parameters.kappa, &overlap_parameters.gamma);

    printf("dt ?\n");
    fgets(renglon, sizeof(renglon), data_file);
    sscanf(renglon, "%f", &overlap_parameters.dt);

    printf("time_run ?\n");
    fgets(renglon, sizeof(renglon), data_file);
    sscanf(renglon, "%f ", &time_run);

    printf("idum ?\n");
    fgets(renglon, sizeof(renglon), data_file);
    sscanf(renglon, "%d", &idum);
}

void print_parameters_eliminate_overlap(FILE *fp, OverlapParameters overlap_parameters, const float time_run){

    fprintf(fp, "%e %e kappa, gamma\n", overlap_parameters.kappa, overlap_parameters.gamma);
    fprintf(fp, "%e dt\n", overlap_parameters.dt);
    fprintf(fp, "%f time_run\n", time_run);
    fprintf(fp, "%d, idum\n", overlap_parameters.idum);
}

void scan_parameters_simulation(FILE *fp, SystemParameters &pars, char *infile, float &trans_time, float &run_time, int &nsamples) {

    char renglon[200];

    printf("temp_set ?\n");
    fgets(renglon, sizeof(renglon), fp);
    sscanf(renglon, "%f", &pars.temperature);

    printf("dt ?\n");
    fgets(renglon, sizeof(renglon), fp);
    sscanf(renglon, "%f", &pars.dt);

    printf("trans_time, run_time, nsamples (taken in run)\n");
    fgets(renglon, sizeof(renglon), fp);
    sscanf(renglon, "%f %f %d", &trans_time, &run_time, &nsamples);

    printf("infile ?\n");
    fgets(renglon, sizeof(renglon), fp);
    sscanf(renglon, "%s", infile);

    printf("nbins_gder, range_gder ?\n");
    fgets(renglon, sizeof(renglon), fp);
    sscanf(renglon, "%d %f", &pars.nbins_gder, &pars.range_gder);

    printf("nbins_dplt, range_dplt ?\n");
    fgets(renglon, sizeof(renglon), fp);
    sscanf(renglon, "%d %f", &pars.nbins_dplt, &pars.range_dplt);
}

void print_parameters_simulation(FILE *fp, const SystemParameters pars, const char *infile, const float trans_time, const float run_time, const int nsamples){

    fprintf(fp, "%f  temperature\n", pars.temperature);
    fprintf(fp, "%f  dt\n", pars.dt);
    fprintf(fp, "%f  %f  %d  trans_time, run_time, nsamples (taken in run)\n", trans_time, run_time, nsamples);
    fprintf(fp, "%s  infile\n", infile);
    fprintf(fp, "%d %f  nbins_gder, range_gder\n", pars.nbins_gder, pars.range_gder);
    fprintf(fp, "%d %f  nbins_dplt, range_dplt\n", pars.nbins_dplt, pars.range_dplt);
}