# NMRRelaxation.jl

This repository has the codes to compute the NMR dipole-dipole autocorrelations for bulk fluids, ideas that formed the basis of the references noted below. The following examples are included for reference: (1) bulk water modeled using the TIP4P potential, and (2) bulk n-heptane described by the CGenFF potential. The water simulations are based on equilibration in the NVT ensemble followed by a production in the NVE ensemble. The n-heptane simulations used NVT ensemble simulations, with the thermostatting performed using the less-intrusive stochastic velocity rescaling thermostat. 


The original codes that informed the references noted below were based on Python, with the parts that needed to be accelerated written in Fortran. The present code is written entirely in Julia. The Julia language approach here is very much that of an eager amateur, so more will need to be done to have it in a way that better conforms to the standard in other Julia packages.  

It should be quite easy to extend the ideas presented here to fluids next to interfaces, as was already done in our published work for heptane in kerogens. 

## Getting started

The NMR relaxation code depends on the [MDToolBox.jl code](https://github.com/matsunagalab/MDToolbox.jl) with couple of changes. The CUDA calls within the MDToolBox code were disabled (since I do not have access to a CUDA device and the only functionality that I needed for the NMR code was being able to extract the trajectory).  I also edited the ``fileIO.jl`` file to return the box dimensions expicitly. Specifically, if we did
```
ta = mdload("w256.psf")
ta = mdload("w256.dcd",top=ta)
```
then I could not obtain the box dimensions which should, in principle, be obtained as
```
box_dimensions = ta.boxsize;
```

I made a simple fix to the ``mdload`` function within  ``fileIO.jl``. Specifically, 

```
   # Code fragment I added 
   if sizeof(ta.boxsize) == 0
       println("No box information available")
       flag = 0
    else
       flag = 1
       boxsize = ta.boxsize
    end

    # Original code 
    if !isnothing(top)
        if isnothing(index)
            ta = [top[0, :]; ta]
        else
            ta = [top[0, index]; ta]
        end
    end

    # Code fragment I added/edited
    if flag == 0
       return ta
    else
       return ta, boxsize
    end
```
This revised package, I call ``MDTools``. 

## Examples 

I have provided two examples. 
1) A 256 molecule box of TIP4P water molecules equilibrated using the CSVR thermostat for 100,000 steps. Then the production phase is in the NVE ensemble for 819200 steps with a frequency of 100 fs to save the frames. Thus we have 8192 = 2<sup>13</sup> frames. (The power of 2 ensures that FFT for convolutions is fast.) The time-step for integration was 1 fs. The temperature of the system was 293.15 K. The data we get matches the results in [Singer, Asthagiri, Chapman, Hirasaki, Journal of Magnetic Resonance, 2017](10.1016/j.jmr.2017.02.001). 
2)  The second example is of n-heptane. The system comprises 261 molecules on n-heptane, and the conditions are exactly those described in the JMR paper noted above. As in the JMR paper I simulated the system for 2 ns (with 1 fs timestep). The one change from the JMR paper was that now I used the velocity rescaling thermostat (with a 1 ps coupling). The results are again consistent with earlier published data. 

## References

If you use the codes or are inspired by it, please cite:

1) P. M. Singer, D. Asthagiri, W. G. Chapman, and G. J. Hirasaki, “Molecular dynamics simulations of
NMR relaxation and diffusion of bulk hydrocarbons and water”, [J. Magn. Reson. 277, 15–24 (2017)](https://dx.doi.org/10.1016/j.jmr.2017.02.001). 
2) P. M. Singer, D. Asthagiri, Z. Chen, A. Valiya Parambathu, G. J. Hirasaki, and W. G. Chapman,
“Role of internal motions and molecular geometry on the NMR relaxation of hydrocarbons”, 
[J. Chem. Phys. 148, 164507 (2018)](https://dx.doi.org/10.1063/1.5023240).
3) D. Asthagiri, W. G. Chapman, G. J. Hirasaki, and P. M. Singer, “NMR 1H-1H dipole relaxation in
fluids: Relaxation of individual 1H-1H pairs versus relaxation of molecular modes”, 
[J. Phys. Chem. B 124, 10802–10810 (2020)](https://dx.doi.org/10.1021/acs.jpcb.0c08078).

## To do

Refactoring the [spin-rotation relaxation](https://dx.doi.org/10.1063/1.5027097) calculation code to Julia. 
