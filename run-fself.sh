#!/bin/bash

#########################Parameters#########################

#set_init_config values

distribution="discrete"
number_species=5
ngrain=1024

phi=(0.20  0.10  0.02  0.01  0.0025)
diameter=(1.00  0.50  0.33  0.25  0.20)
polydispersity=(0.00  0.10  0.10  0.05  0.05)
max_polydispersity=(0.00  0.50  0.50  0.25  0.25)

#eliminateOverlap values
kappa=50000
gamma=200
dt_overlap=5.0e-4
runtime_overlap=4.0

#creating configs values
dt_configs=1.0e-4
runtime_configs=0.1
transtime_configs=0.1
number_configs=5
tagb=8
tags=8

#Temperature values
type="linear"
temp0=1.0
tempf=0.25
t0=0.0
tf=0.0
ngap_rescaling=1

#simulation parameters
sampling="log"      # sampling type: linear or log
dt_aging=1.0e-4
transtime_aging=0.1
runtime_aging=1.0e1
initial_decade=1.0e-3   # only for log sampling
number_initial_decades=3
samples_aging=10


#g(r) and f_self(r)

gr_nbins=400
gr_range=15.0
qmax=6.19

# verify parameters for set_init_config

if [ "${#phi[@]}" -ne "$number_species" ]; then
  echo "Error, number of phi values does not match number of species"
  exit
fi

if [ "${#diameter[@]}" -ne "$number_species" ]; then
  echo "Error, number of diameter values does not match number of species"
  exit
fi

if [ "${#polydispersity[@]}" -ne "$number_species" ]; then
  echo "Error, number of polydispersity values does not match number of species"
  exit
fi

if [ "${#max_polydispersity[@]}" -ne "$number_species" ]; then
  echo "Error, number of max_polydispersity values does not match number of species"
  exit
fi

#Crea carpeta para el sistema particular

if [ $# -eq 0 ]; then
	echo "Error, no se ingresaron argumentos"
	exit 1
fi

directorio=ns_${number_species}_${distribution}_temp_${temp0}_${tempf}_${tf}_sample$1

if [ -d "$directorio" ]; then
	echo "Directorio ya existe, use otro numero"
	exit 2
fi

mkdir $directorio
cp run-fself.sh $directorio/
cd $directorio



###################Crea carpetas ###################
echo "Creating files"

for((i=0; i<=number_initial_decades;i++)); do
  mkdir -p "configs/decade$i"
done

for((i=0; i<=number_initial_decades;i++)); do
  mkdir -p "simulation/results/decade$i"
done

cd configs

#copy files to simulation
cp ../../setConfig.bin .
cp ../../eliminateOverlap.bin .
cp ../../bd_aging_configs.bin .

#create setConfig.data
idum1=$RANDOM
{
echo "$distribution"
echo "$number_species"

for phi_value in "${phi[@]}"; do
  echo -n "$phi_value  "
done
echo

for diameter_value in "${diameter[@]}"; do
  echo -n "$diameter_value  "
done
echo

for polydispersity_value in "${polydispersity[@]}"; do
  echo -n "$polydispersity_value  "
done
echo

for max_polydispersity_value in "${max_polydispersity[@]}"; do
  echo -n "$max_polydispersity_value  "
done
echo

echo "$ngrain"
echo "1.0"     #temperature
echo "$idum1"
echo "init_config_overlaped"
} > setConfig.data

# set values in eliminate_overlap.data
{
  echo "init_config_overlaped"
  echo "rcp_bitacora    rcp_energias"
  echo "rcp_snapshots   init_config"
  echo "$kappa  $gamma"
  echo "$dt_overlap"
  echo "$runtime_overlap"
  echo "1.0 $RANDOM"
} > eliminateOverlap.data

# set values in wca-configs.data
{
  echo "$temp0"
  echo "$dt_configs"
  echo "$transtime_configs  $runtime_configs  $number_configs"
  echo "$tagb  $tags"
  echo "init_config"
} > bd_aging_configs.data

cd ..

#copy files to simulation
cp ../bd_aging_self.bin simulation/

###########Init configs############

cd configs

#creating init_config
echo "Creating initial config"
./setConfig.bin < setConfig.data #> bitacora_set.log
./eliminateOverlap.bin < eliminateOverlap.data #> bitacora_overlap.log
#rm init_config_overlaped
#rm rcp*

#creating all starting configs
echo "Creating all configs"
./bd_aging_configs.bin < bd_aging_configs.data #> bitacora_configs.log

mv init_config_* decade0

cd ..

###########simulation############
echo "Starting simulation"

cd simulation

for(( i=0; i<=number_initial_decades; i++ )); do

  if [ $i -eq 0 ]; then
    print_decades_flag=1
  else
    print_decades_flag=0
  fi
  
  current_decade=$i

  {
    echo "$type"

    if [ "$type" = "linear" ]; then
      echo "$temp0  $tempf  $t0  $tf  $ngap_rescaling"
    elif [ "$type" = "sine" ]; then
      exit 3 # not implemented yet
    fi

    echo "$dt_aging"
    echo "$sampling"

    if [ "$sampling" = "linear" ]; then
        echo "$transtime_aging  $runtime_aging  $samples_aging"
      elif [ "$sampling" = "log" ]; then
        echo "$transtime_aging  $initial_decade  $runtime_aging"
      else
        echo "sampling option not implemented"
        exit 3 # not implemented yet
    fi

    echo "$print_decades_flag  $number_initial_decades  $current_decade"
    echo "$tagb  $tags"
    echo "$number_configs"
    echo "$gr_nbins  $gr_range"
    echo "$qmax"
    echo "$RANDOM"
  } > bd_aging_self.data


  #run simulation
  ./bd_aging_self.bin < bd_aging_self.data #> bitacora_runs.log

  if [ $? -ne 0 ]; then
  echo "Error in simulation"
  exit 4
  fi

done

echo "Simulation finished successfully"
