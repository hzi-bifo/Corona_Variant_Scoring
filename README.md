### Corona Variant Prediction

Pipeline to score pangolin lineages based on their number of known antigenic sites. 

Will require the following repositories to be installed:

Corona_Antigenic_Sites

To run the pipeline you will need to specify the input/output directory as the first argument and change the path to the Corona_Antigenic_Sites pipeline and software path in the variant_scoring.sh script. The input directory most contain a GISAID metadata .tsv file that will be used as the initila input. Please see the test/ directory for an example.

Run the following command to start the pipeline:

bash variant_scoring.sh <"path to desired input / output directory">
