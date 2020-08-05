# Tuolumne basin image mode/classic mode comparison

Files are in `/home/jschap/Documents/Codes/VICMATLAB/examples/`. Demo the subsetting codes and publish the instructions in a Markdown file, which is available for reference from the VICGlobal Zenodo respository. Compare Classic/Image mode outputs for Tuolumne basin.

Simulation uses the Livneh et al. (2013) meteorological forcing data, and is for the calendar years 2009-2011, with no accounting for soil moisture spin-up. Energy balance mode with five snowbands. 

See MATLAB wrapper for the comparison. Will write this up here as a Markdown file, once it is running. 

<!-- Run the classic mode simulation -->
```bash
./vic_classic.exe -g /home/jschap/Documents/Codes/VICMATLAB/data/VICGlobal/global_param.txt
```

<!-- Run the image mode simulation -->
```bash
./vic_classic.exe -g /home/jschap/Documents/Codes/VICMATLAB/data/VICGlobal/global_param.txt
```