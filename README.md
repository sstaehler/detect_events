# detect_events
Detect (synthetic) marsquakes in real noise from Mars.


# Installation

## Prequesites
The game makes use of [ObsPy](https://www.obspy.org) and [Instaseis](https://www.instaseis.net). Please make sure to have these packages installed. 

    conda install -c conda-forge matplotlib obspy
    pip install instaseis

## Download and install

    git clone https://github.com/sstaehler/detect_events
    cd detect_events
    pip install -v -e .
    

## Noise time series
The code is delivered with a four day file of synthetic martian noise from the MQS Operational Readiness Test. 
If you prefer to use real Martian noise, you need to be a team member. If you are part of the InSight project, contact 
Simon Stähler for the password to download it from [here](https://polybox.ethz.ch/index.php/s/jvs2DbB1MKTlmlV) and save
it to 


    detect_events/data/
    
    




